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

#define MIN_VALUE 0.000001

using namespace geos::geom;
using namespace geos::algorithm;

class line {
public:
	double x0;
	double y0;
	double x1;
	double y1;
	line(double x0, double y0, double x1, double y1) :
		x0(x0), y0(y0), x1(x1), y1(y1) {}
};

class inpoint {
public:
	double w; // Wavelength
	double ss; // Sample spectra
	inpoint(double w, double ss) :
		w(w), ss(ss) {}
};

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



class input {
public:
	int c, r;
	std::vector<inpoint> data;
	input() : input(0, 0) {}
	input(int c, int r) :
		c(c), r(r) {
	}
};

class output {
public:
	int c, r;
	double area; // Hull area.
	double larea; // Left hand area
	double rarea; // Right hand area.
	double symmetry; // larea / area
	std::vector<outpoint> data;
	output() : output(0, 0) {}
	output(input& in) :
		output(in.c, in.r) {}
	output(int c, int r) :
		c(c), r(r),
		area(0), larea(0), rarea(0), symmetry(0) {}
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

	const GeometryFactory gf;

	// Make a list of Coordinates.
	std::vector<Coordinate> coords;
	for(const inpoint& pt : in)
		coords.emplace_back(pt.w, pt.ss, 0);

	// Make a MultiPoint from the coords and get the ConvexHull.
	MultiPoint* mp = gf.createMultiPoint(coords);
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

void processQueue(std::list<input>* inqueue, std::list<output>* outqueue,
		std::mutex* inmtx, std::condition_variable* incv, std::mutex* outmtx, std::condition_variable* outcv, std::condition_variable* readcv, bool* running) {

	input in;
	while(true) {
		{
			std::unique_lock<std::mutex> lk(*inmtx);
			// If the input is empty or the output is too large, pause.
			while((inqueue->empty() || outqueue->size() > 1000) && *running)
				incv->wait(lk);
			// If input is empty and reading is done, quit the loop.
			if(inqueue->empty() && !*running)
				break;
			in = inqueue->front();
			inqueue->pop_front();
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
		double maxCrm = 0, maxCr = 0;
		int maxIdx = 0, i = 0;
		for(outpoint& pt : out.data) {
			pt.cr = pt.ss / pt.ch;
			pt.crm = 1 - pt.cr;
			if(pt.cr > maxCr)
				maxCr = pt.cr;
			if(pt.crm > maxCrm) {
				maxCrm = pt.crm;
				maxIdx = i;
			}
			++i;
		}

		// Calculate the other ch metrics.
		for(outpoint& pt : out.data) {
			pt.crn = pt.cr / maxCr;
			pt.crm = 1 - pt.cr;
			pt.crnm = pt.crm / maxCrm;
		}

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
		out.rarea = out.area - out.larea;
		out.symmetry = out.larea / out.rarea;

		{
			// Send to output queue and notify.
			std::lock_guard<std::mutex> lk(*outmtx);
			outqueue->push_back(out);
		}

		outcv->notify_one();
		if(inqueue->size() < 1000)
			readcv->notify_one();

		//std::cerr << "inqueue " << inqueue->size() << "; outqueue " << outqueue->size() << "\n";

	}
}

void writeQueue(const std::string* outfile, const std::vector<std::string>& wavelengths, int cols, int rows, int bands, std::list<output>* outqueue,
		std::mutex* outmtx, std::condition_variable* outcv, std::condition_variable* incv, bool* running) {

	GDALWriter writerss(*outfile + "_ss.tif", cols, rows, bands, "wavelength", wavelengths);
	GDALWriter writerch(*outfile + "_ch.tif", cols, rows, bands, "wavelength", wavelengths);
	GDALWriter writercr(*outfile + "_cr.tif", cols, rows, bands, "wavelength", wavelengths);
	GDALWriter writercrn(*outfile + "_crn.tif", cols, rows, bands, "wavelength", wavelengths);
	GDALWriter writercrm(*outfile + "_crm.tif", cols, rows, bands, "wavelength", wavelengths);
	GDALWriter writercrnm(*outfile + "_crnm.tif", cols, rows, bands, "wavelength", wavelengths);
	GDALWriter writerhull(*outfile + "_hull.tif", cols, rows, 4, "stat", {"hull_area", "hull_left_area", "hull_right_area", "hull_symmetry"});

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
			std::unique_lock<std::mutex> lk(*outmtx);
			while(outqueue->empty() && *running)
				outcv->wait(lk);
			if(outqueue->empty() && !*running)
				break;
			out = outqueue->front();
			outqueue->pop_front();
		}

		for(const outpoint& o : out.data) {
			ss.push_back(o.ss);
			ch.push_back(o.ch);
			cr.push_back(o.cr);
			crn.push_back(o.crn);
			crm.push_back(o.crm);
			crnm.push_back(o.crnm);
		}

		hull = {out.area, out.larea, out.rarea, out.symmetry};

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

		incv->notify_one();

		//std::cerr << "outqueue " << outqueue->size() << "\n";
	}

	writerhull.writeStats(*outfile + "_hullstats.csv", {"hull_area", "hull_left_area", "hull_right_area", "hull_symmetry"});

}

void Processor::process(Reader& reader, const std::string& outfile, int bufSize, int threads) {

	std::list<input> inqueue;
	std::list<output> outqueue;
	std::mutex inmtx;
	std::mutex outmtx;
	std::mutex readmtx;
	std::condition_variable incv;
	std::condition_variable outcv;
	std::condition_variable readcv;
	bool inrunning = true;
	bool outrunning = true;

	std::vector<double> buf(bufSize * bufSize * reader.bands());

	const std::vector<double>& wavelengths = reader.getBands();

	std::vector<std::string> wavelengthMeta;
	for(double w : wavelengths)
		wavelengthMeta.push_back(std::to_string(w));

	std::list<std::thread> t0;
	for(int i = 0; i < threads; ++i)
		t0.emplace_back(processQueue, &inqueue, &outqueue, &inmtx, &incv, &outmtx, &outcv, &readcv, &inrunning);
	std::thread t1(writeQueue, &outfile, wavelengthMeta, reader.cols(), reader.rows(), reader.bands(), &outqueue, &outmtx, &outcv, &incv, &outrunning);


	int col, row, cols, rows;
	while(reader.next(buf, col, row, cols, rows)) {

		std::cerr << col << "," << row << " of " << reader.cols() << "," << reader.rows() << "\n";

		{
			// If input queue gets too large, wait before adding more.
			std::unique_lock<std::mutex> lk(readmtx);
			if(inqueue.size() > 1000)
				readcv.wait(lk);
		}

		{
			// Read out the input objects and add to the queue.
			std::lock_guard<std::mutex> lk(inmtx);
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c) {
					input in(c + col, r + row);
					for(int b = 0; b < reader.bands(); ++b) {
						double v = buf[b * bufSize * bufSize + r * bufSize + c];
						double w = wavelengths[b];
						in.data.emplace_back(w, v);
					}
					inqueue.push_back(in);
				}
			}
		}

		// Notify the input processor.
		incv.notify_all();
	}

	inrunning = false;
	incv.notify_all();

	for(std::thread& t : t0)
		t.join();

	outrunning = false;
	outcv.notify_all();

	t1.join();
}


