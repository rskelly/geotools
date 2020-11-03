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

	// Handle the arguments.
	std::string outfile = argv[1];
	std::vector<std::string> infiles;
	for(int i = 2; i < argc; ++i)
		infiles.push_back(argv[i]);

	// Get the number of bands in the first file.
	int bandCount = Band<float>::bands(infiles[0]);

	geo::util::Bounds<double> bounds;
	GridProps props;
	std::vector<std::string> bandMetaNames;
	std::vector<std::vector<std::string>> bandMetadata;

	g_debug("Getting band descriptors...")
	for(int i = 0; i < bandCount; ++i) {
		Band<float> band(infiles.front(), i, false, false);
		bandMetaNames.push_back(band.props().bandMetaName());
		bandMetadata.push_back(band.props().bandMetadata());
	}

	g_debug("Calculating bounds...");
	for(size_t i = 0; i < infiles.size(); ++i) {
		const std::string& file = infiles[i];
		int b;
		if((b = Band<float>::bands(file)) != bandCount)
			g_runerr("Files  must all have the same number of bands. " << bandCount << " expected; " << b << " found in " << file);
		Band<float> band(file, 0, false, false);
		bounds.extend(band.props().bounds());
		if(i == 0)
			props = band.props();
	}

	// Extend the props object using the overall bands and make writable.
	props.setBounds(bounds);
	props.setWritable(true);

	// To keep track of counts.
	std::vector<int> counts(props.cols() * props.rows());
	std::fill(counts.begin(), counts.end(), 0);

	g_debug("Making temporary files...");
	std::vector<Band<float>*> bands;
	for(int i = 0; i < bandCount; ++i) {
		GridProps p(props);
		p.setBandMetaName(bandMetaNames[i]);
		p.setBandMetadata(bandMetadata[i]);	
		bands.push_back(new Band<float>("/tmp/band_" + std::to_string(i) + ".tif", p, true));
		bands[i]->fill(0);
	}

	g_debug("Adding values...");
	for(const std::string& file : infiles) {

		// Get the count band.
		Band<float> cband(file, 0, false, false);
		const GridProps& p = cband.props();
		int cols = p.cols();
		int rows = p.rows();
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				int cc = props.toCol(p.toX(c));
				int rr = props.toRow(p.toY(r));
				// Sum the counts.
				if(cband.get(c, r) != props.nodata())
					counts[rr * props.cols() + cc] += cband.get(c, r);
			}
		}
		
		for(int i = 1; i < bandCount; ++i) {
			Band<float> band(file, i, false, false);
			for(int r = 0; r < rows; ++r) {
				for(int c = 0; c < cols; ++c) {
					int cc = props.toCol(p.toX(c));
					int rr = props.toRow(p.toY(r));
					// Add the weighted value to the main band value.
					if(band.get(c, r) != props.nodata())
						bands[i]->set(cc, rr, bands[i]->get(cc, rr) + counts[rr * props.cols() + cc] * band.get(c, r));
				}
			}
		}
	}

	g_debug("Normalizing...");
	int cols = props.cols();
	int rows = props.rows();
	for(int i = 1; i < bandCount; ++i) {
		Band<float>* band = bands[i];
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				// Normalize against the count.
				if(counts[r * cols + c]) {
					band->set(c, r, band->get(c, r) / counts[r * cols + c]);
				} else {
					band->set(c, r, props.nodata());
				}
			}
		}
	}
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {
			// Set the zero counts to nodata.
			bands[0]->set(c, r, counts[r * cols + c] ? counts[r * cols + c] : props.nodata());
		}
	}

	Band<float>::mergeBands(bands, outfile, "GTiff", true);

}


