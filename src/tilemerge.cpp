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
	Band<int> counts(props, true);
	Band<float> sums(props, true);
	counts.fill(0);

	g_debug("Making temporary files...");
	std::vector<Band<float>*> bands;
	for(int i = 0; i < bandCount; ++i) {
		GridProps p(props);
		p.setBandMetaName(bandMetaNames[i]);
		p.setBandMetadata(bandMetadata[i]);	
		bands.push_back(new Band<float>("/tmp/band_" + std::to_string(i) + ".tif", p, true));
		bands[i]->fill(0);
	}


	g_debug("Adding counts.");
	for(const std::string& file : infiles) {
		Band<float> cband(file, 0, false, false);
		const GridProps& p = cband.props();
		for(int r = 0; r < p.rows(); ++r) {
			for(int c = 0; c < p.cols(); ++c) {
				float v = cband.get(c, r);
				if(v > 0 && v != p.nodata()) {
					int cc = props.toCol(p.toX(c));
					int rr = props.toRow(p.toY(r));
					counts.set(cc, rr, counts.get(cc, rr) + (int) v);
				}
			}
		}
	}

	g_debug("Calculating.");
	for(int i = 1; i < bandCount; ++i) {
		sums.fill(0);
		for(const std::string& file : infiles) {
			Band<float> cband(file, 0, false, false);
			Band<float> vband(file, i, false, false);
			const GridProps& p = vband.props();
			for(int r = 0; r < p.rows(); ++r) {
				for(int c = 0; c < p.cols(); ++c) {
					float v = vband.get(c, r);
					int ct = cband.get(c, r);
					if(ct && v != p.nodata()) {
						int cc = props.toCol(p.toX(c));
						int rr = props.toRow(p.toY(r));
						sums.set(cc, rr, sums.get(cc, rr) + v * ct);
					}
				}
			}
		}
		for(int r = 0; r < props.rows(); ++r) {
			for(int c = 0; c < props.cols(); ++c) {
				int ct = counts.get(c, r);
				if(ct) {
					bands[i]->set(c, r, sums.get(c, r) / ct);
				} else {
					bands[i]->set(c, r, props.nodata());
				}
			}
		}
	}
	for(int r = 0; r < props.rows(); ++r) {
		for(int c = 0; c < props.cols(); ++c) {
			int ct = counts.get(c, r);
			if(ct) {
				bands[0]->set(c, r, ct);
			} else {
				bands[0]->set(c, r, props.nodata());
			}
		}
	}

	Band<float>::mergeBands(bands, outfile, "GTiff", true);

}


