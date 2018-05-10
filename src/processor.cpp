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

#include <geos/algorithm/ConvexHull.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/MultiPoint.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Point.h>

#include "processor.hpp"

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
	double crm; // Mirrored cr
	double crnm; // Mirrored normalized cr
	outpoint(inpoint& in, double ch) :
		w(in.w), ss(in.ss),
		ch(ch), cr(0), crm(0), crnm(0) {}
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
	if(x < x0 && x > x1) {
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

	const GeometryFactory* gf = GeometryFactory::getDefaultInstance();

	// Make a list of Coordinates.
	std::vector<Coordinate> coords;
	for(const inpoint& pt : in)
		coords.emplace_back(pt.w, std::max(MIN_VALUE, pt.ss));

	// Make a MultiPoint from the coords and get the ConvexHull.
	MultiPoint* mp = gf->createMultiPoint(coords);
	Polygon* hull = dynamic_cast<Polygon*>(mp->convexHull());

	// Get the area for output.
	area = hull->getArea();

	// Extract the line segments.
	std::vector<line> lines;
	const LineString* ring = hull->getExteriorRing();
	for(size_t i = 0; i < ring->getNumPoints() - 1; ++i) {
		Point* p0 = ring->getPointN(i);
		Point* p1 = ring->getPointN(i + 1);
		double x0 = p0->getX(), y0 = p0->getY();
		double x1 = p1->getX(), y1 = p1->getY();
		// Only add a segment if it isn't at a corner where y=0.
		if(y0 > 0.0 && y1 > 0.0)
			lines.emplace_back(x0, y0, x1, y1);
	}

	return lines;
}

void processQueue(std::list<input>* queue, std::mutex* mtx) {

	input in;
	while(true) {
		{
			std::lock_guard<std::mutex> lk(*mtx);
			in = queue->front();
			queue->pop_front();
		}

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

		double maxCrm = 0;
		int maxIdx = 0, i = 0;
		for(outpoint& pt : out.data) {
			pt.cr = pt.ss / pt.ch;
			pt.crm = 1 - pt.cr;
			if(pt.crm > maxCrm) {
				maxCrm = 0;
				maxIdx = i;
			}
			++i;
		}
		for(outpoint& pt : out.data) {
			pt.crm = 1 - pt.cr;
			pt.crnm = pt.crm / maxCrm;
		}

		// Add two corner points to make the hull full.
		pts.assign(in.data.begin(), in.data.begin() + maxIdx);
		pts.emplace_back(pts[pts.size() - 1].w, 0.0);
		pts.emplace_back(pts[0].w, 0.0);

		// Compute the hull, assign area to the output.
		lines = convexHull(pts, out.larea);

		// Compute the rest of the numbers.
		out.rarea = out.area - out.larea;
		out.symmetry = out.larea / out.area;

	}
}

void Processor::process(const std::vector<double>& buffer, const std::vector<double>& wavelengths,
		std::vector<double>& output, int cols, int rows, int bands, int bufSize) {

	if(wavelengths.size() != (size_t) bands)
		std::invalid_argument("Wavelengths array must be the same size as the number of bands.");

	std::list<input> queue;
	std::mutex mtx;

	std::thread t(processQueue, &queue, &mtx);

	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {
			input in(c, r);
			for(int b = 0; b < bands; ++b) {
				double v = buffer[b * bufSize * bufSize + r * cols + c];
				double w = wavelengths[0];
				in.data.emplace_back(w, v);
			}
			queue.push_back(in);
		}
	}

	t.join();

}


