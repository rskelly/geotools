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

#include "matplotlibcpp.h"

#include "contrem.hpp"
#include "reader.hpp"
#include "writer.hpp"
#include "contrem_util.hpp"

using namespace geos::geom;
using namespace geos::algorithm;
using namespace hlrg;

namespace {

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
		double w; 		// Wavelength.
		double ss; 		// Sample spectra (intensity).
		double ch; 		// Intersection with convex hull (y).
		double cr; 		// Continuum removal (ss/ch).
		double crn;		// Continuum removal normalized (1 - cr).
		double crm;		// Mirrored cr.
		double crnm;	// Mirrored normalized cr.

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
		double maxWl; // Wavelength at the maxCrm value
		double slope;	// The slope of the regression line through the crnm points.
		double yint;	// The y-intercept of the regression.
		int maxCount; // The number of equal maximum values.
		std::vector<outpoint> data;
		output() : output(0, 0) {}
		output(input& in) :
			output(in.c, in.r) {}
		output(int c, int r) :
			c(c), r(r),
			area(0), larea(0), rarea(0), symmetry(0),
			maxCrm(0), maxWl(0),
			slope(0), yint(0),
			maxCount(0) {}
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
			// Only add a segment if it isn't a bottom or end segment (which have at least one zero y-coordinate).
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
		Contrem* contrem;
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

		// Check for a mask file.
		bool hasRoi = false;
		std::unique_ptr<GDALReader> mask;
		{
			std::string roi = config->contrem->roi;
			if(!roi.empty()) {
				hasRoi = true;
				mask.reset(new GDALReader(roi));
			}
		}

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

			// If there's a mask, check it. Skip if necessary.
			if(hasRoi && mask->getInt(in.c, in.r) <= 0)
				continue;

			// Adjust <=0 intensities to MIN_VALUE. This enables the creation
			// of a hull even though the area of the hull will be zero for practical purposes.
			for(size_t i = 0; i < in.data.size(); ++i) {
				if(in.data[i].ss <= MIN_VALUE)
					in.data[i].ss = MIN_VALUE;
			}

			// Add two corner points to complete the hull.
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
				pt.cr = pt.ss / pt.ch;		// Continuum removal -- intensity proportional to height of hull.
				pt.crm = 1 - pt.cr;			// Mirrored continuum removal.
				if(pt.crm > out.maxCrm) {
					// There may be more than one equal max; it'll be flagged at a later step
					out.maxCrm = pt.crm;
					out.maxWl = pt.w;
					maxIdx = i;
				}
				++i;
			}

			double sxy = 0, sx = 0, sy = 0, sxx = 0;
			size_t n = out.data.size();

			// Calculate the other ch metrics.
			for(size_t i = 0; i < n; ++i) {
				outpoint& pt = out.data[i];
				pt.crn = pt.crm / out.maxCrm;	// Normalized continuum removal -- normalized against maximum mirrored cr.
				pt.crnm = 1 - pt.crn;			// Mirrored normalized continuum removal.
				if(i > 0 && i < n - 1) {
					sxy += pt.w * pt.crnm;
					sx += pt.w;
					sy += pt.crnm;
					sxx += pt.w * pt.w;
				}
				// Count the values equal to max; will flag those with >1.
				if(pt.crm == out.maxCrm)
					++out.maxCount;
			}

			out.slope = ((n - 2) * sxy - sx * sy) / ((n - 2) * sxx - sx * sx);
			out.yint = (sy - out.slope * sx) / (n - 2);

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

		std::string outfile = config->contrem->output;
		std::string driver = config->contrem->outputType;
		const std::vector<double>& wavelengths = config->wavelengths;
		const std::vector<std::string>& bandNames = config->bandNames;
		int cols = config->cols;
		int rows = config->rows;
		int bands = config->bands;

		std::string ext;
		if(driver == "ENVI") {
			ext = "";
		} else if(driver == "GTiff") {
			ext = ".tif";
		} else {
			throw std::invalid_argument("Unknown driver: " + driver);
		}

		// Remove the extension if there is one.
		outfile = outfile.substr(0, outfile.find_last_of("."));

		char* meta;

		GDALWriter writerss(outfile + "_ss" + ext, driver, cols, rows, bands, wavelengths, bandNames);
		GDALWriter writerch(outfile + "_ch" + ext, driver, cols, rows, bands, wavelengths, bandNames);
		GDALWriter writercr(outfile + "_cr" + ext, driver, cols, rows, bands, wavelengths, bandNames);
		GDALWriter writercrnm(outfile + "_crnm" + ext, driver, cols, rows, bands, wavelengths, bandNames);
		GDALWriter writerhull(outfile + "_agg" + ext, driver, cols, rows, 9, {}, {"hull_area", "hull_left_area", "hull_right_area", "hull_symmetry", "max_crm", "max_crm_wl", "max_count", "slope", "y-int"});
		GDALWriter writermax(outfile + "_maxcount" + ext, driver, cols, rows, 1, {}, {"equal_max_count"}, &meta, DataType::Byte);
		GDALWriter writervalid(outfile + "_valid" + ext, driver, cols, rows, 1, {}, {"valid_hull"}, &meta, DataType::Byte);

		// Each row is a string of coords which can be turned into a hull.
		std::ofstream hulls(outfile + "_hulls.csv");
		hulls << "col,row,slope,yint\n";

		writermax.fill(0);
		writervalid.fill(0);

		std::vector<double> ss;
		std::vector<double> ch;
		std::vector<double> cr;
		std::vector<double> crnm;
		std::vector<double> hull;
		std::vector<double> w;

		std::vector<int> maxima(cols * rows);
		std::vector<int> valid(cols * rows);

		std::fill(maxima.begin(), maxima.end(), 0);
		std::fill(valid.begin(), valid.end(), 0);

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

			hulls << out.c << "," << out.r << "," << out.slope << "," << out.yint << "\n";

			for(const outpoint& o : out.data) {
				ss.push_back(o.ss);
				ch.push_back(o.ch);
				cr.push_back(o.cr);
				crnm.push_back(o.crnm);
				w.push_back(o.w);
			}
			hulls << "\n";

			{
				namespace plt = matplotlibcpp;
				std::vector<double> rx = {w.front(), w.back()};
				std::vector<double> ry = {w.front() * out.slope + out.yint, w.back() * out.slope + out.yint};
				plt::named_plot("Normalized, Continuum Removed", w, crnm);
				plt::named_plot("Regression", rx, ry);
				plt::title(std::string("Normalized Continuum Removal, ") + std::to_string(out.c) + "," + std::to_string(out.r));
				plt::save(outfile + "/hull_img/hull_" + std::to_string(out.c) + "_" + std::to_string(out.r) + ".png");
				plt::close();
			}

			maxima[out.r * cols + out.c] = out.maxCount > 1 ? 0 : 1;
			valid[out.r * cols + out.c] = out.area > 0 && out.rarea > 0 && out.larea > 0;

			hull = {out.area, out.larea, out.rarea, out.symmetry, out.maxCrm, out.maxWl, (double) out.maxCount, out.slope, out.yint};

			writerss.write(ss, out.c, out.r, 1, 1, 1, 1);
			writerch.write(ch, out.c, out.r, 1, 1, 1, 1);
			writercr.write(cr, out.c, out.r, 1, 1, 1, 1);
			writercrnm.write(crnm, out.c, out.r, 1, 1, 1, 1);
			writerhull.write(hull, out.c, out.r, 1, 1, 1, 1);

			ss.clear();
			ch.clear();
			cr.clear();
			crnm.clear();

			config->incv.notify_one();
		}

		writermax.write(maxima, 0, 0, cols, rows);
		writervalid.write(valid, 0, 0, cols, rows);

		writerhull.writeStats(outfile + "_agg_stats.csv", {"hull_area", "hull_left_area", "hull_right_area", "hull_symmetry", "max_crm", "max_crm_wl", "max_count", "slope", "yint"});
	}

}

std::unique_ptr<Reader> getReader(const std::string& file, bool transpose, int headerRows, int minCol, int maxCol) {
	std::unique_ptr<Reader> rdr;
	FileType type = getFileType(file);
	switch(type) {
	case CSV:
		rdr.reset(new CSVReader(file, transpose, headerRows, minCol, maxCol));
		break;
	case GTiff:
	case ENVI:
		rdr.reset(new GDALReader(file));
		break;
	case ROI:
	case SHP:
	case SQLITE:
	case Unknown:
	default:
		throw std::runtime_error("Unknown file type for " + file);
	}
	return rdr;
}

void Contrem::run(ContremListener* listener) {

	QConfig qconfig;
	qconfig.contrem = this;

	std::unique_ptr<Reader> reader = getReader(spectra, wlTranspose, wlHeaderRows, wlMinCol, wlMaxCol);

	reader->setBandRange(minWl, maxWl);

	qconfig.cols = reader->cols();
	qconfig.rows = reader->rows();
	qconfig.bands = reader->bands();

	// A buffer for input data. Stores a row from a raster, or a single "pixel" from a table.
	std::vector<double> buf(qconfig.cols * reader->bands());

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
	for(int i = 0; i < threads; ++i)
		t0.emplace_back(processQueue, &qconfig);

	// Start the output thread.
	std::thread t1(writeQueue, &qconfig);

	// Read through the buffer and populate the input queue.
	int cols, row;
	int bands = reader->bands();
	while(reader->next(buf, cols, row)) {

		{
			// Read out the input objects and add to the queue.
			std::lock_guard<std::mutex> lk(qconfig.inmtx);
			for(int c = 0; c < cols; ++c) {
				input in(c, row);
				for(int b = 0; b < bands; ++b) {
					double v = buf[c * bands + b];
					double w = qconfig.wavelengths[b];
					in.data.emplace_back(w, v);
				}
				qconfig.inqueue.push_back(in);
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

double Contrem::progress() const {
	return 0;
}

/*
 * 		if(m_reader)
			delete m_reader;
		if(!m_roiFile.empty()) {
			m_reader = new ROIReader(m_roiFile);
		} else if(!m_spectraFile.empty()) {
			m_reader = new GDALReader(m_spectraFile);
		} else {
			throw std::invalid_argument("No input file (-r or -d) given.");
		}

		if(m_minWl > 0 && m_maxWl > 0)
			m_reader->setBandRange(m_minWl, m_maxWl);

		m_reader->setBufSize(m_buffer);
 *
 */
