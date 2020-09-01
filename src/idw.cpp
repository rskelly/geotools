/*
 * splinesmooth.cpp
 *
 *  Created on: Oct 27, 2019
 *      Author: rob
 */

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

#include "util.hpp"
#include "grid.hpp"

using namespace geo::util;
using namespace geo::util::csv;
using namespace geo::grid;

class Pt {
public:
	double _x, _y, _z;

	Pt(double x = 0, double y = 0, double z = 0) :
		_x(x), _y(y), _z(z) {}

	void x(double x) {
		_x = x;
	}
	void y(double y) {
		_y = y;
	}
	void z(double z) {
		_z = z;
	}
	double x() const {
		return _x;
	}
	double y() const {
		return _y;
	}
	double z() const {
		return _z;
	}
	double& operator[](int idx) {
		switch(idx) {
		case 0: return _x;
		case 1: return _y;
		default: return _z;
		}
	}
};
void usage() {
	std::cout << "Usage: idw [options] <points> <columns> <outfile>\n"
			" -rx <res x>        The output grid resolution in x.\n"
			" -ry <res y>        The output grid resolution in y.\n"
			" -s <srid>          The projection of the output grid.\n"
			" -b <buffer>        A buffer around the maxima of the point set to define\n"
			"                    the bounds of the output raster.\n"
			" -t <raster>        A template raster. Supercedes the resolution, projection,"
			"                    srid and buffer parameters.\n"
			" -h                 If there's a header in the csv point file, use this switch.\n\n"
			" -m <smooth>        The smoothing parameter. If not given or less than or equal to zero, \n"
			"                    the number of input points is used.\n"
			" -d <denominator>   Decimate the points. The given number is the denominator of the retained fraction.\n"
			" <points>           Is a CSV file containing at least x, y and z columns with zero or one header lines.\n"
			" <columns>          A comma-delimited list of indices of columns in the csv file\n"
			"                    for the x, y and z columns. An optional fourth column will be\n"
			"                    used for weights. This must be accompanied by the -m switch with\n"
			"                    a smoothing factor. If weights are not given, the std. deviation\n"
			"                    of the z-coordinates is used.\n"
			" <outfile>          The name of a geotiff to write output to.\n";
}

int main(int argc, char** argv) {

	std::vector<std::string> args;

	double xres = 0, yres = 0;		// Output resolution if no template is given.
	int srid = 0;					// Spatial reference ID for the output (must match csv file).
	double buffer = 0;				// Buffer to add to the point bounds if no template is given.
	std::string tpl;				// Filename of template raster.
	std::string projection;			// Raster output projection.
	bool header = false;			// True if there is one field header in the csv file.
	std::vector<int> columns;		// Column indices for the csv file.
	double minx = geo::maxvalue<double>(), miny = geo::maxvalue<double>(),
			maxx = geo::minvalue<double>(), maxy = geo::minvalue<double>();	// The dataset bounds. If there's a template, use those bounds, otherwise use the buffered point bounds.
	bool hasTemplate = false;		// Set to true if a template is loaded.
	double smooth = 0;				// The smoothing parameter for the bivariate spline.
	bool csv = false;
	double decimate = 1;

	for(int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if(arg == "-rx") {
			xres = atof(argv[++i]);
		} else if(arg == "-ry") {
			yres = atof(argv[++i]);
		} else if(arg == "-s") {
			srid = atoi(argv[++i]);
		} else if(arg == "-b") {
			buffer = atof(argv[++i]);
		} else if(arg == "-t") {
			tpl = argv[++i];
		} else if(arg == "-h") {
			header = true;
		} else if(arg == "-m") {
			smooth = atof(argv[++i]);
			if(smooth < 0)
				smooth = 0;
		} else if(arg == "-d") {
			decimate = 1. / atof(argv[++i]);
		} else {
			args.push_back(argv[i]);
		}
	}

	if(args.size() < 3) {
		std::cerr << "Too few arguments.\n";
		usage();
		return 1;
	}

	if(args[2].find(".csv") != std::string::npos) {
		std::cout << "CSV filename given. Using CSV mode.\n";
		csv = true;
	}

	if(!csv && tpl.empty() && (xres == 0 || yres == 0)) {
		std::cerr << "If a template is not given, and not in CSV mode, xres and yres must be nonzero.\n";
		usage();
		return 1;
	}

	// Parse the column indices.
	{
		std::vector<std::string> cs;
		split(std::back_inserter(cs), args[1], ",");
		for(const std::string& c : cs)
			columns.push_back(atoi(c.c_str()));
		if(columns.size() < 3) {
			std::cerr << "Too few csv columns. There must be three.\n";
			usage();
			return 1;
		}
	}

	GridProps props;				// Properties for the output grid.

	// Load the template raster.
	if(!tpl.empty()) {
		try {
			Band<float> tplg(tpl, 0, false, true);
			const GridProps tprops = tplg.props();
			const Bounds<double>& tbounds = tprops.bounds();
			xres = tprops.resX();
			yres = tprops.resY();
			projection = tprops.projection();
			minx = tbounds.minx();
			miny = tbounds.miny();
			maxx = tbounds.maxx();
			maxy = tbounds.maxy();
			props = tprops;
			hasTemplate = true;
		} catch(const std::runtime_error& ex) {
			std::cerr << "Failed to load template raster.\n";
			usage();
			return 1;
		}
	}

	// Load the CSV data.
	std::vector<Pt> pts;
	{
		CSV csv(args[0], header);
		std::vector<CSVValue> cx = csv.column(columns[0]);
		std::vector<CSVValue> cy = csv.column(columns[1]);
		std::vector<CSVValue> cz = csv.column(columns[2]);
		std::vector<CSVValue> cw;
		if(columns.size() > 3)
			cw = csv.column(columns[3]);
		if(!(cx.size() == cy.size() && cy.size() == cz.size())) {
			std::cerr << "Input coordinate arrays must be the same length.\n";
			usage();
			return 1;
			if(!cw.empty() && cx.size() != cw.size()) {
				std::cerr << "The weights list must be the same length as the coordinate arrays.\n";
				usage();
				return 1;
			}
		}
		for(size_t i = 0; i < cx.size(); ++i)
			pts.emplace_back(cx[i].asDouble(), cy[i].asDouble(), cz[i].asDouble());
	}

	if(decimate < 1.) {
		std::cout << "Decimating: " << decimate << "\n";
		std::vector<size_t> idx(pts.size());
		for(size_t i = 0; i < pts.size(); ++i)
			idx[i] = i;
		std::random_shuffle(idx.begin(), idx.end());
		size_t s = std::max(1, (int) std::ceil(pts.size() * decimate));
		std::vector<Pt> _pts(s);
		for(size_t i = 0; i < s; ++i)
			_pts[i] = std::move(pts[idx[i]]);
		pts.swap(_pts);
	}

	// Set the bounds if there's no template.
	if(!hasTemplate) {
		minx = miny = std::numeric_limits<double>::max();
		maxx = maxy = std::numeric_limits<double>::lowest();

		for(const Pt& p : pts) {
			if(p.x() < minx) minx = p.x();
			if(p.x() > maxx) maxx = p.x();
			if(p.y() < miny) miny = p.y();
			if(p.y() > maxy) maxy = p.y();
		}

		minx -= buffer;
		miny -= buffer;
		maxx += buffer;
		maxy += buffer;

	}

	if(csv) {

		std::ofstream csv(args[2]);
		csv << "x,y,z\n";

		for(size_t i = 0; i < pts.size(); ++i) {
			double x = pts[i].x();
			double y = pts[i].z();
			double s = 0;
			double w = 0;
			for(size_t j = 0; j < pts.size(); ++j) {
				double d = 1.0 / std::pow(pts[j].x() - x, 2.0) + std::pow(pts[j].y() - y, 2.0);
				s += pts[j].z() * d;
				w += d;
			}
			csv << pts[i].x() << "," << pts[i].y() << "," << (s / w) << "\n";
		}

	} else {
		std::cout << "Starting raster...\n";

		props.setResolution(xres, yres);
		props.setProjection(projection);
		props.setSrid(srid);
		props.setBounds(Bounds<double>(minx, miny, maxx, maxy));
		props.setNoData(-9999.0);
		props.setDataType(DataType::Float32);
		props.setWritable(true);
		props.setBands(1);

		Band<float> outgrid(args[2], props);

		int cols = props.cols();
		int rows = props.rows();

		double sigma = 500.0;
		double norm = 1.0 / (sigma * std::sqrt(2 * M_PI));

		for(int r = 0; r < rows; ++r) {
			std::cout << "Row " << r << " of " << rows << "\n";
			for(int c = 0; c < cols; ++c) {
				double x = props.toX(c);
				double y = props.toY(r);
				double s = 0;
				double w = 0;
				for(size_t j = 0; j < pts.size(); ++j) {
					double d0 = std::sqrt(std::pow(pts[j].x() - x, 2.0) + std::pow(pts[j].y() - y, 2.0));
					double w0 = norm * std::exp(-0.5 * std::pow(d0 / sigma, 2.0));
					s += pts[j].z() * w0;
					w += w0;
				}
				if(w != 0) {
					outgrid.set(x, y, s / w);
				} else {
					outgrid.set(x, y, props.nodata());
				}
			}
		}
	}


}

