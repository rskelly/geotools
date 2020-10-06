/*
 * srtmlakes.cpp
 *
 * Locates flat regions in SRTM DEMs and polygonizes them, under
 * the assumption that they're lakes. A minimum area threshold prevents
 * the inclusion of small clusters of equal-value pixels.
 *
 *  Created on: Sep 21, 2020
 *      Author: rob
 */

#include <vector>
#include <string>

#include "geo.hpp"
#include "util.hpp"
#include "grid.hpp"

using namespace geo::util;
using namespace geo::grid;

void usage() {
	std::cout << "Usage: srtmlakes <min area> <multiplier> <sqlite file> <in file [in file [in file ...]]>\n";
}

int main(int argc, char** argv) {

	int minarea = 0;
	float multiplier = 1.0;
	std::string outfile;
	std::vector<std::string> infiles;

	for(int i = 1; i < argc; ++i) {
		if(i == 1) {
			minarea = atof(argv[i]);
		} else if(i == 2){
			multiplier = atof(argv[i]);
		} else if(i == 3) {
			outfile = argv[i];
		} else if(i > 3) {
			infiles.emplace_back(argv[i]);
		}
	}

	if(outfile.empty()) {
		g_error("No output file given.");
		usage();
		return 1;
	}

	if(infiles.empty()) {
		g_error("No input files given.");
		usage();
		return 1;
	}

	if(minarea <= 0) {
		g_error("Minimum area must be larger than zero.");
		usage();
		return 1;
	}

	if(multiplier <= 0) {
		g_error("Multiplier must be larger than zero.");
		usage();
		return 1;
	}

	GridProps props;
	Bounds<double> bounds;
	Band<int> input;
	int nodata;
	{
		for(size_t i = 0; i < infiles.size(); ++i) {
			Band<int16_t> grd(infiles[i], 0, false, true);
			bounds.extend(grd.props().bounds());
			if(i == 0)
				props = grd.props();
		}
		props.setBounds(bounds);
		props.setWritable(true);
		props.setBands(1);
		props.setDataType(DataType::Int32);
		nodata = props.nodata();
		input.init(geo::util::tmpfile("srtmlakes"), props, true);
		for(size_t i = 0; i < infiles.size(); ++i) {
			Band<int16_t> grd(infiles[i], 0, false, true);
			const GridProps& pr = grd.props();
			grd.writeTo(input, pr.cols(), pr.rows(),
					0, 0,
					props.toCol(pr.toX(0)), props.toRow(pr.toY(0))
			);
		}
		if(multiplier != 1) {
			float v;
			for(int r = 0; r < props.rows(); ++r) {
				for(int c = 0; c < props.cols(); ++c) {
					if((v = input.get(c, r)) != nodata) {
						input.set(c, r, (int) std::round(v * multiplier));
					} else {
						input.set(c, r, nodata);
					}
				}
			}
		}
	}

	GridProps lakeProps(props);
	lakeProps.setNoData(0);
	Band<int> lakes(geo::util::tmpfile("srtmlakes"), lakeProps, true);
	lakes.fill(0);

	TargetFillOperator<int, int> op1(&input, &lakes, 0, -1);			// Fills from target to -1.
	TargetFillOperator<int, int> op2(&lakes, &lakes, -1, 0);			// Fills from -1 to final value.
	TargetFillOperator<int, int> op3(&input, &input, 0, nodata);		// Fills from target to nodata.

	int minc, minr, maxc, maxr, area;
	int vi;
	int id = 0;
	float cellArea = std::abs(props.resX() * props.resY());

	int statusCount = std::max(1, props.rows() / 10);
	for(int r = 0; r < props.rows(); ++r) {
		if(r % statusCount == 0)
			std::cout << "Row " << r << " of " << props.rows() << "\n";
		for(int c = 0; c < props.cols(); ++c) {
			if((vi = input.get(c, r)) != nodata && vi >= 0) {
				op1.setTarget(vi);
				Band<int>::floodFill(c, r, op1, false, &minc, &minr, &maxc, &maxr, &area);
				if(area * cellArea >= minarea) {
					op2.setFill(++id);
					Band<int>::floodFill(c, r, op2, false, &minc, &minr, &maxc, &maxr, &area);
				}
				op3.setTarget(vi);
				Band<int>::floodFill(c, r, op3, false, &minc, &minr, &maxc, &maxr, &area);
			}
		}
	}

	lakes.polygonizeToFile(outfile, "lakes", "lakeid", "SQLite", props.projection());

}


