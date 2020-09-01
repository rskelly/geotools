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

void usage() {
	std::cout << "Usage: splinesmooth [options] <points> <columns> <outfile>\n"
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
			" -f                 Force overwriting of existing files.\n"
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
	bool force = false;

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
		} else if(arg == "-f") {
			force = true;
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

	if(!tpl.empty() && !isfile(tpl)) {
		g_error("A template file is given but it does not exist.");
		usage();
		return 1;
	}

	if(!geo::util::safeToWrite(args[2], force)) {
		g_error("The file " << args[2] << " is not safe to write to. "
				<< "If it is a file, use the -f flag. If it is a directory, "
				<< "choose a different path.");
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
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	std::vector<double> w;
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
		for(const CSVValue& v : cx)
			x.push_back(v.asDouble());
		for(const CSVValue& v : cy)
			y.push_back(v.asDouble());
		for(const CSVValue& v : cz)
			z.push_back(v.asDouble());
		if(!cw.empty()) {
			for(const CSVValue& v : cw)
				w.push_back(v.asDouble());
		}

	}

	if(decimate < 1.) {
		std::cout << "Decimating: " << decimate << "\n";
		std::vector<size_t> idx(x.size());
		for(size_t i = 0; i < x.size(); ++i)
			idx[i] = i;
		std::random_shuffle(idx.begin(), idx.end());
		size_t s = std::max(1, (int) std::ceil(x.size() * decimate));
		std::vector<double> _x(s), _y(s), _z(s), _w;
		bool dow = w.size() == x.size();
		if(dow)
			_w.resize(s);
		for(size_t i = 0; i < s; ++i) {
			_x[i] = x[idx[i]];
			_y[i] = y[idx[i]];
			_z[i] = z[idx[i]];
			if(dow)
				_w[i] = w[idx[i]];
		}
		_x.swap(x);
		_y.swap(y);
		_z.swap(z);
		if(dow)
			_w.swap(w);
	}

	// Set the bounds if there's no template.
	if(!hasTemplate) {
		minx = miny = std::numeric_limits<double>::max();
		maxx = maxy = std::numeric_limits<double>::lowest();

		for(double xx : x) {
			if(xx < minx) minx = xx;
			if(xx > maxx) maxx = xx;
		}
		for(double yy : y) {
			if(yy < miny) miny = yy;
			if(yy > maxy) maxy = yy;
		}

		minx -= buffer;
		miny -= buffer;
		maxx += buffer;
		maxy += buffer;

	}

	std::cout << "Preparing bivariate spline with " << x.size() << " coordinates.\n";
	BivariateSpline bvs;
	bvs.init(smooth, x, y, z, w, minx, miny, maxx, maxy);

	w.clear();

	if(csv) {

		std::ofstream csv(args[2]);
		csv << "x,y,z\n";

		std::vector<double> xx(1);
		std::vector<double> yy(1);
		std::vector<double> zz(1);
		for(size_t i = 0; i < x.size(); ++i) {
			xx[0] = x[i];
			yy[0] = y[i];
			bvs.evaluate(xx, yy, zz);
			csv << xx[0] << "," << yy[0] << "," << zz[0] << "\n";
		}

	} else {
		std::cout << "Starting raster...\n";

		props.setResolution(xres, yres);
		props.setProjection(projection);
		props.setSrid(srid);
		props.setBounds(Bounds(minx, miny, maxx, maxy));
		props.setNoData(-9999.0);
		props.setDataType(DataType::Float32);
		props.setWritable(true);
		props.setBands(1);

		Band<float> outgrid(args[2], props);

		int cols = props.cols();
		int rows = props.rows();
		x.resize(1);
		y.resize(1);
		z.resize(1);
		for(int r = 0; r < rows; ++r) {
			std::cout << "Row " << r << " of " << rows << "\n";
			for(int c = 0; c < cols; ++c) {
				x[0] = props.toX(c);	// TODO: This only seems to work one cell at a time...
				y[0] = props.toY(r);
				bvs.evaluate(x, y, z);
				outgrid.set(x[0], y[0], z[0]);
			}
		}
	}


}

