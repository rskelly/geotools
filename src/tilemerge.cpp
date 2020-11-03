/*
 * tilemerge.cpp
 *
 * This program merges tiled output from pc2grid into a single raster.
 *
 * It uses the point count on the first layer to perform weighted
 * averaging of values on subsequent layers to eliminate edge effects
 * where the tiles overlap.
 *
 *  Created on: Nov 3, 2020
 *      Author: rob
 */

#include <string>

#include "geo.hpp"
#include "grid.hpp"
#include "util.hpp"

using namespace geo::grid;
using namespace geo::util;

int main(int argc, char** argv) {

	if(argc < 3) {
		g_error("Usage: tilemerge <outfile> <infile [infile [...]]>");
		return 1;
	}

	std::string outfile = argv[1];
	std::vector<std::string> infiles;
	for(int i = 2; i < argc; ++i)
		infiles.push_back(argv[i]);

	g_debug("Calculating bounds...");
	geo::util::Bounds<double> bounds;
	GridProps props;
	for(size_t i = 0; i < infiles.size(); ++i) {
		const std::string& file = infiles[i];
		Band<float> band(file, 0, false, false);
		bounds.extend(band.props().bounds());
		if(i == 0)
			props = band.props();
	}

	props.setBounds(bounds);
	props.setWritable(true);

	std::vector<int> counts(props.cols() * props.rows());
	std::fill(counts.begin(), counts.end(), 0);

	g_debug("Making temporary files...");
	std::vector<Band<float>*> bands;
	for(int i = 0; i < props.bands(); ++i)
		bands.push_back(new Band<float>("/tmp/band_" + std::to_string(i) + ".tif", props, true));

	g_debug("Adding values...");
	for(const std::string& file : infiles) {
		Band<float> cband(file, 0, false, false);
		for(int i = 0; i < props.bands(); ++i) {
			Band<float> band(file, i, false, false);
			const GridProps& p = band.props();
			int cols = p.cols();
			int rows = p.rows();
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c) {
					double x = p.toX(c);
					double y = p.toY(r);
					int cc = props.toCol(x);
					int rr = props.toRow(y);
					// Add the weighted value to the main band value.
					bands[i]->set(cc, rr, bands[i]->get(cc, rr) + cband.get(c, r) * band.get(c, r));
					counts[rr * props.cols() + cc] += cband.get(c, r);
				}
			}
		}
	}

	g_debug("Normalizing...");
	for(Band<float>* band : bands) {
		int cols = props.cols();
		int rows = props.rows();
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				// Normalize against the count.
				band->set(c, r, band->get(c, r) / counts[r * cols + c]);
			}
		}
	}

	Band<float>::mergeBands(bands, outfile, "GTiff", true);

}


