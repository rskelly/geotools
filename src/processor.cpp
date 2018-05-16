/*
 * processor.cpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#include <vector>
#include <list>
#include <iostream>
#include <mutex>
#include <unordered_map>
#include <cmath>
#include <thread>
#include <condition_variable>

#include <geos/algorithm/ConvexHull.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/MultiPoint.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Point.h>

#include "processor.hpp"
#include "reader.hpp"
#include "writer.hpp"

using namespace geos::geom;
using namespace geos::algorithm;

/**
 * A line, representing a segment from a convex hull.
 */
class line {
public:
	double x0;
	double y0;
	double x1;
	double y1;
	line(double x0, double y0, double x1, double y1) :
		x0(x0), y0(y0), x1(x1), y1(y1) {}
};

/**
 * An input point. Contains a wavelength and sample intensity.
 */
class inpoint {
public:
	double w; // Wavelength
	double ss; // Sample spectra
	inpoint(double w, double ss) :
		w(w), ss(ss) {}
};

/**
 * An output point. Contains the same info as an input point,
 * plus continuum removal values.
 */
class outpoint {
public:
	double w; // Wavelength
	double ss; // Sample spectra
	double ch; // Intersection with convex hull (y)
	double cr; // Continuum removal (ss/ch)
	double crn;
	double crm; // Mirrored cr
	double crnm; // Mirrored normalized cr
	outpoint(inpoint& in, double ch) :
		w(in.w), ss(in.ss),
		ch(ch), cr(0), crn(0), crm(0), crnm(0) {}
	outpoint(inpoint& in) :
		outpoint(in, 0) {}
};

/**
 * An input pixel, which includes spectral information
 * and cell coordinates.
 */
class input {
public:
	int c, r;
	std::vector<inpoint> data;
	input() : input(0, 0) {}
	input(int c, int r) :
		c(c), r(r) {
	}
};

/**
 * An output pixel which contains convex hull information,
 * continuum removal information and cell coordinates.
 */
class output {
public:
	int c, r;
	double area; // Hull area.
	double larea; // Left hand area
	double rarea; // Right hand area.
	double symmetry; // larea / area
	double maxCrm; // maximum normalized mirrored continuum removal
	std::vector<outpoint> data;
	output() : output(0, 0) {}
	output(input& in) :
		output(in.c, in.r) {}
	output(int c, int r) :
		c(c), r(r),
		area(0), larea(0), rarea(0), symmetry(0), maxCrm(0) {}
};

/**
 * Returns the y value corresponding to x along the given line.
 * NaN if no intersection found.
 */
double interpolate(double x, double x0, double y0, double x1, double y1) {
	if(x < x0 || x > x1) {
		return std::numeric_limits<double>::quiet_NaN();
	} else if(x == x0) {
		return y0;
	} else if(x == x1) {
		return y1;
	} else {
		return y0 + (x - x0) / (x1 - x0) * (y1 - y0);
	}
}

/**
 * Compute the convex hull around the points and return the line segments.
 */
std::vector<line> convexHull(const std::vector<inpoint>& in, double& area) {

	const GeometryFactory::unique_ptr gf = GeometryFactory::create();

	// Make a list of Coordinates.
	std::vector<Coordinate> coords;
	for(const inpoint& pt : in)
		coords.emplace_back(pt.w, pt.ss, 0);

	// Make a MultiPoint from the coords and get the ConvexHull.
	MultiPoint* mp = gf->createMultiPoint(coords);
	Polygon* hull = dynamic_cast<Polygon*>(mp->convexHull());

	// Get the area for output.
	area = hull->getArea();

	// Extract the line segments.
	std::vector<line> lines;
	const LineString* ring = hull->getExteriorRing();
	for(size_t i = 0; i < ring->getNumPoints() - 1; ++i) {
		Point* p0 = ring->getPointN(i); // This call implies an allocation. LAME. Delete at the end.
		Point* p1 = ring->getPointN(i + 1);
		double x0 = p0->getX(), y0 = p0->getY();
		double x1 = p1->getX(), y1 = p1->getY();
		// Only add a segment if it isn't at a corner where y=0.
		if(y0 > 0 && y1 > 0)
			lines.emplace_back(x0, y0, x1, y1);
		delete p0;
		delete p1;
	}

	delete mp;
	delete hull;

	return lines;
}

class QConfig {
public:
	ProcessorConfig pconfig;
	std::list<input> inqueue;
	std::list<output> outqueue;
	std::mutex inmtx;
	std::mutex outmtx;
	std::mutex readmtx;
	std::condition_variable incv;
	std::condition_variable outcv;
	std::condition_variable readcv;
	bool inRunning;
	bool outRunning;

	int cols;
	int rows;
	int bands;
	std::vector<double> wavelengths;
	std::vector<std::string> bandNames;

	std::unordered_map<size_t, bool> maximaFlag;

	std::vector<std::string> getWavelengthNames() const {
		std::vector<std::string> names;
		for(double w : wavelengths)
			names.push_back(std::to_string(w));
		return names;
	}
};

/**
 * Process the input queue.
 */
void processQueue(QConfig* config) {

	input in;
	while(true) {
		{
			std::unique_lock<std::mutex> lk(config->inmtx);
			// If the input is empty or the output is too large, pause.
			while((config->inqueue.empty() || config->outqueue.size() > 1000) && config->inRunning)
				config->incv.wait(lk);
			// If input is empty and reading is done, quit the loop.
			if(config->inqueue.empty() && !config->inRunning)
				break;
			in = config->inqueue.front();
			config->inqueue.pop_front();
		}

		// Adjust <=0 intensities to MIN_VALUE.
		for(size_t i = 0; i < in.data.size(); ++i)
			if(in.data[i].ss <= MIN_VALUE) in.data[i].ss = MIN_VALUE;

		// Add two corner points to make the hull full.
		std::vector<inpoint> pts(in.data.begin(), in.data.end());
		pts.emplace_back(pts[pts.size() - 1].w, 0.0);
		pts.emplace_back(pts[0].w, 0.0);

		// Create an output pixel to hold computed values.
		output out(in);

		// Compute the hull, assign area to the output.
		std::vector<line> lines = convexHull(pts, out.area);

		// Find the intersection point for each wavelength.
		for(inpoint& pt : in.data) {
			bool found = false;
			for(line& l : lines) {
				double ch = interpolate(pt.w, l.x0, l.y0, l.x1, l.y1);
				if(!std::isnan(ch)) {
					out.data.emplace_back(pt, ch);
					found = true;
					break;
				}
			}
			if(!found)
				throw std::runtime_error("Failed to find an intersection point.");
		}

		// Calculate the cr and crm, and get the max value and index.
		out.maxCrm = 0;
		int maxIdx = 0, i = 0;
		for(outpoint& pt : out.data) {
			pt.cr = pt.ss / pt.ch;
			pt.crm = 1 - pt.cr;
			if(pt.crm > out.maxCrm) {
				out.maxCrm = pt.crm;
				maxIdx = i;
			}
			++i;
		}

		// Calculate the other ch metrics.
		int maxCount = 0;
		for(size_t i = 0; i < out.data.size(); ++i) {
			outpoint& pt = out.data[i];
			pt.crn = pt.crm / out.maxCrm;
			pt.crnm = 1 - pt.crn;
			if(pt.crm > out.maxCrm)
				out.maxCrm = pt.crnm;
			if(pt.crm == out.maxCrm)
				++maxCount;
		}

		if(maxCount > 1) {
			// There is more than one equal maximum. Flag it.
			config->maximaFlag[((size_t) out.c << 16) | out.r] = true;
		}

		// We were going to do interpolation for adjacent maxima, but put it off.
		// Kopăcková, V., & Koucká, L. (2017). Integration of absorption feature information from visible to
		// longwave infrared spectral ranges for mineral mapping. Remote Sensing, 9(10), 8–13. https://doi.org/10.3390/rs9101006
		// If there are  >2 maxima, or the distance between them is > than the configured
		// interp distance, flag the cell and move on. Otherwise, interpolate.

		// Calculate the split hull; add two corner points to make the hull full.
		pts.assign(in.data.begin(), in.data.begin() + maxIdx);

		if(pts.size() > 1) {
			pts.emplace_back(pts[pts.size() - 1].w, 0.0);
			pts.emplace_back(pts[0].w, 0.0);

			// Compute the hull, assign area to the left output.
			lines = convexHull(pts, out.larea);
		} else {
			out.larea = 0;
		}

		// Compute the rest of the numbers.
		if(out.area == 0 || out.larea == 0 || out.larea == out.area) {
			// If the left or right area is zero, or the total is zero, we
			out.larea = 0;
			out.area = 0;
			out.symmetry = 0;
			out.rarea = 0;
		} else {
			out.rarea = out.area - out.larea;
			out.symmetry = out.larea / out.rarea;
		}

		{
			// Send to output queue and notify.
			std::lock_guard<std::mutex> lk(config->outmtx);
			config->outqueue.push_back(out);
		}

		config->outcv.notify_one();
		config->readcv.notify_one();

	}
}

/**
 * Process the output queue and write to file.
 */
void writeQueue(QConfig* config) {

	const std::string& outfile = config->pconfig.outfile;
	const std::string& driver = config->pconfig.driver;
	const std::string& ext = config->pconfig.extension;
	const std::vector<double>& wavelengths = config->wavelengths;
	const std::vector<std::string>& bandNames = config->bandNames;
	int cols = config->cols;
	int rows = config->rows;
	int bands = config->bands;

	GDALWriter writerss(outfile + "_ss" + ext, driver, cols, rows, bands, wavelengths, bandNames);
	GDALWriter writerch(outfile + "_ch" + ext, driver, cols, rows, bands, wavelengths, bandNames);
	GDALWriter writercr(outfile + "_cr" + ext, driver, cols, rows, bands, wavelengths, bandNames);
	GDALWriter writercrn(outfile + "_crn" + ext, driver, cols, rows, bands, wavelengths, bandNames);
	GDALWriter writercrm(outfile + "_crm" + ext, driver, cols, rows, bands, wavelengths, bandNames);
	GDALWriter writercrnm(outfile + "_crnm" + ext, driver, cols, rows, bands, wavelengths, bandNames);
	GDALWriter writerhull(outfile + "_hull" + ext, driver, cols, rows, 5, {}, {"hull_area", "hull_left_area", "hull_right_area", "hull_symmetry", "max_crnm"});
	GDALWriter writermax(outfile + "_maxima" + ext, driver, cols, rows, 1, {}, {"maximum"}, DataType::Byte);

	std::vector<double> ss;
	std::vector<double> ch;
	std::vector<double> cr;
	std::vector<double> crn;
	std::vector<double> crm;
	std::vector<double> crnm;
	std::vector<double> hull;

	output out;
	while(true) {
		{
			std::unique_lock<std::mutex> lk(config->outmtx);
			while(config->outqueue.empty() && config->outRunning)
				config->outcv.wait(lk);
			if(config->outqueue.empty() && !config->outRunning)
				break;
			out = config->outqueue.front();
			config->outqueue.pop_front();
		}

		for(const outpoint& o : out.data) {
			ss.push_back(o.ss);
			ch.push_back(o.ch);
			cr.push_back(o.cr);
			crn.push_back(o.crn);
			crm.push_back(o.crm);
			crnm.push_back(o.crnm);
		}

		hull = {out.area, out.larea, out.rarea, out.symmetry, out.maxCrm};

		writerss.write(ss, out.c, out.r, 1, 1, 1);
		writerch.write(ch, out.c, out.r, 1, 1, 1);
		writercr.write(cr, out.c, out.r, 1, 1, 1);
		writercrn.write(crn, out.c, out.r, 1, 1, 1);
		writercrm.write(crm, out.c, out.r, 1, 1, 1);
		writercrnm.write(crnm, out.c, out.r, 1, 1, 1);
		writerhull.write(hull, out.c, out.r, 1, 1, 1);

		ss.clear();
		ch.clear();
		cr.clear();
		crn.clear();
		crm.clear();
		crnm.clear();

		config->incv.notify_one();
	}

	writerhull.writeStats(outfile + "_stats.csv", {"hull_area", "hull_left_area", "hull_right_area", "hull_symmetry", "max_crm"});

	std::vector<int> maxima(cols * rows);
	for(const auto& item : config->maximaFlag) {
		int c = (item.first >> 16) & 0xffff;
		int r = item.first & 0xffff;
		maxima[r * cols + c] = item.second;
	}
	writermax.write(maxima, 0, 0, cols, rows, 1);
}

void Processor::process(Reader* reader, const ProcessorConfig& config) {

	QConfig qconfig;
	qconfig.pconfig = config;
	qconfig.cols = reader->cols();
	qconfig.rows = reader->rows();
	qconfig.bands = reader->bands();

	int bufSize = config.bufferSize;

	// A buffer for input data.
	std::vector<double> buf(bufSize * bufSize * reader->bands());

	// A list of wavelengths.
	qconfig.wavelengths = reader->getWavelengths();
	qconfig.bandNames = reader->getBandNames();

	qconfig.inRunning = true;
	qconfig.outRunning = true;

	// A list of wavelengths as strings for labelling.
	std::vector<std::string> wavelengthMeta;
	for(double w : qconfig.wavelengths)
		wavelengthMeta.push_back(std::to_string(w));

	// Start the processing threads.
	std::list<std::thread> t0;
	for(int i = 0; i < config.threads; ++i)
		t0.emplace_back(processQueue, &qconfig);

	// Start the output thread.
	std::thread t1(writeQueue, &qconfig);

	// Read through the buffer and populate the input queue.
	int col, row, cols, rows;
	while(reader->next(buf, col, row, cols, rows)) {

		std::cerr << col << "," << row << " of " << reader->cols() << "," << reader->rows() << "\n";

		{
			// Read out the input objects and add to the queue.
			std::lock_guard<std::mutex> lk(qconfig.inmtx);
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c) {
					input in(c + col, r + row);
					for(int b = 0; b < reader->bands(); ++b) {
						double v = buf[b * bufSize * bufSize + r * bufSize + c];
						double w = qconfig.wavelengths[b];
						in.data.emplace_back(w, v);
					}
					qconfig.inqueue.push_back(in);
				}
			}
		}

		// Notify the input processor.
		qconfig.incv.notify_all();

		// If input queue gets too large, wait before adding more.
		std::unique_lock<std::mutex> lk(qconfig.readmtx);
		while(qconfig.inqueue.size() > 1000)
			qconfig.readcv.wait(lk);
	}

	// Let the processor threads finish.
	qconfig.inRunning = false;
	qconfig.incv.notify_all();
	for(std::thread& t : t0)
		t.join();

	// Let the output thread finish.
	qconfig.outRunning = false;
	qconfig.outcv.notify_all();
	t1.join();
}


