/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include <fstream>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <string>
#include <vector>

#include "geo.hpp"
#include "pointcloud.hpp"
#include "util.hpp"
#include "grid.hpp"

using namespace geo;
using namespace geo::pc;
using namespace geo::util;
using namespace geo::grid;

void usage() {
	std::cerr << "Usage: pc2grid [options] <output raster> <input las [*]>\n"
			<< "If the output raster is to be merged (default), the given filename \n"
			<< " is used. Otherwise, each output band is saved to a file named "
			<< "<basename>_<method>.<extension>.\n"
			<< "\n"
			<< " -r <filename>    Use a template raster for resolution, extent and projection.\n"
			<< "                  Overrides rx/ry, e, n and s.\n"
			<< " -rx <resolution> The output x resolution in map units.\n"
			<< " -ry <resolution> The output y resolution in map units.\n"
			<< " -e <easting>     The top left corner horizontal alignment\n"
			<< "                  (defaults to nearest multiple of resolution).\n"
			<< " -n <northing>    The top left corner vertical alignment\n"
			<< "                  (defaults to nearest whole multiple of resolution).\n"
			<< " -d <radius>      The search radius for finding points. Set to zero\n"
			<< "                  to use the rectangular cell bounds. Defaults to root\n"
			<< "                  of 1/2 cell diagonal (a circle that touches the corners.\n"
			<< " -s <srid>        The spatial reference ID. Default 0.\n"
			<< " -m <methods>     Comma-separated list of statistics to compute; one on each \n"
			<< "                  layer (see below.) Default mean.\n"
			<< "                  rugosity, variance, std. deviation and percentile.\n"
			<< "                  For percentile, use the form, 'percenile:n', where\n"
			<< "                  n is the percentile (no % sign); 1 - 99.\n"
			<< " -v               Verbose. Enable debug and warning messages.\n"
			<< " -h               If given *do not* trust the LAS file headers to contain\n"
			<< "                  good bounds (etc.) info\n"
			<< " -t <min>         Normalize the point density in each cell. Any cell with more than this\n"
			<< "                  number of points will be randomly thinned. Any cell with less will be\n"
			<< "                  reduced to zero. This occurs after filtering.\n"
			<< " -d <nodata>      A nodata value. The default is -9999.0\n"
			<< " -b <bounds>      A comma-delimited list of coordinates, minx, miny, maxx, maxy giving the size of the\n"
			<< "                  raster in projected coordinates.\n"
			<< " -g               Do not merge individual bands to a single file using the given output filename.\n"
			<< " -f               Force the overwrite of existing output files.\n"
			<< " -l               The LRU cache node size. Smaller is faster. Larger uses less memory. Default 1000.\n";

	PCPointFilter::printHelp(std::cerr);

	std::cerr << "\n Available computers: \n";
	for(auto& item : geo::pc::Rasterizer::availableComputers())
		std::cerr << " - " << item.first << ": " << item.second << ".\n";

}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	double resX = std::nan("");
	double resY = std::nan("");
	double easting = std::nan("");
	double northing = std::nan("");
	double radius = std::nan("");
	uint16_t srid = 0;
	bool useHeader = true;
	std::vector<std::string> types;
	std::vector<std::string> args;
	PCPointFilter filter;
	int thin = 0;
	double nodata = -9999;
	double bounds[4] = {std::nan("")};
	std::string templ;
	std::string projection;
	bool merge = true;
	bool force = false;
	int lruSize = 1000;

	for(int i = 1; i < argc; ++i) {
		if(filter.parseArgs(i, argv))
			continue;
		std::string v = argv[i];
		if(v == "-m") {
			std::string type = argv[++i];
			split(std::back_inserter(types), lowercase(type), ",");
		} else if(v == "-r") {
			templ = argv[++i];
		} else if(v == "-h") {
			useHeader = false;
		} else if(v == "-v") {
			geo::loglevel(G_LOG_TRACE);
		} else if(v == "-rx") {
			resX = atof(argv[++i]);
		} else if(v == "-ry") {
			resY = atof(argv[++i]);
		} else if(v == "-d") {
			radius = atof(argv[++i]);
		} else if(v == "-s") {
			srid = atoi(argv[++i]);
		} else if(v == "-e") {
			easting = atof(argv[++i]);
		} else if(v == "-n") {
			northing = atof(argv[++i]);
		} else if(v == "-t") {
			thin = atoi(argv[++i]);
		} else if(v == "-d") {
			nodata = atof(argv[++i]);
		} else if(v == "-f") {
			force = true;
		} else if(v == "-l") {
			lruSize = atoi(argv[++i]);
		} else if(v == "-b") {
			std::vector<std::string> parts;
			split(std::back_inserter(parts), argv[++i]);
			if(parts.size() < 4)
				g_runerr("Not enough parts in the bounds string.");
			for(size_t i = 0; i < 4; ++i) {
				bounds[i] = atof(parts[i].c_str());
			}
		} else if(v == "-g") {
			merge = false;
		} else {
			args.push_back(argv[i]);
		}
	}

	std::string outfile = args.front();
	std::vector<std::string> infiles(args.begin() + 1, args.end());

	// If any files have a * or ? character, run the search.
	{
		std::vector<std::string> tmp;
		for(const std::string& f : infiles) {
			if(f.find("*") < std::string::npos || f.find("?") < std::string::npos) {
				std::vector<std::string> matches = geo::util::glob(f);
				if(matches.empty()) {
					g_error("No files found for search string: " << f);
					usage();
					return 1;
				} else {
					tmp.insert(tmp.end(), matches.begin(), matches.end());
				}
			} else {
				tmp.push_back(f);
			}
		}
		infiles.swap(tmp);
	}

	filter.print();
	
	if(args.size() < 2) {
		g_error("Input and output filenames required.");
		usage();
		return 1;
	}

	if(!templ.empty()) {
		if(!isfile(templ)) {
			g_error("A template file is given but it does not exist.");
			usage();
			return 1;
		}
		g_debug("Getting raster parameters from template: " << templ);
		Band<float> tgrid(templ, 0, false, true);
		const GridProps& props = tgrid.props();
		resX = props.resX();
		resY = props.resY();
		easting = props.tlx();
		northing = props.tly();
		projection = props.projection();
	} else if(srid > 0) {
		projection = projectionFromSRID(srid);
	}

	if(!checkValidInputFiles(infiles)) {
		g_error("At least one of the input files is invalid or doesn't exist.");
		for(const std::string& f : infiles)
			g_error(" - " << f);
		usage();
		return 1;
	}

	if(!geo::util::safeToWrite(outfile, force)) {
		g_error("The file " << outfile << " is not safe to write to. "
				<< "If it is a file, use the -f flag. If it is a directory, "
				<< "choose a different path.");
		usage();
		return 1;
	}

	try {
		Rasterizer r(infiles);
		r.setFilter(&filter);
		r.setThin(thin);
		r.setNoData(nodata);
		r.setBounds(bounds);
		r.setMerge(merge);
		r.setLRUSize(lruSize);
		r.rasterize(outfile, types, resX, resY, easting, northing, radius, projection, useHeader);
	} catch(const std::exception& ex) {
		std::cerr << ex.what() << "\n";
		usage();
		return 1;
	}
	return 0;
}


