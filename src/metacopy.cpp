/*
 * metacopy.cpp
 *
 * Copy band metadata from one raster to another.
 *
 *  Created on: Feb 7, 2020
 *      Author: rob
 */

#include "grid.hpp"

using namespace geo::grid;

int main(int argc, char** argv) {

	if(argc < 5) {
		std::cerr << "Usage: metacopy <source file> <target file> <dest file> <property name>\n";
		return 1;
	}

	std::string srcfile = argv[1];
	std::string tgtfile = argv[2];
	std::string dstfile = argv[3];
	std::string pname = argv[4];

	geo::grid::Grid<float> srcg(srcfile);
	const GridProps& srcp = srcg.props();
	std::vector<std::string> meta = srcp.bandMetadata();

	geo::grid::Grid<float> tgtg(tgtfile);
	const GridProps& tgtp = tgtg.props();

	if(meta.size() > tgtp.bands()) {
		meta.resize(tgtp.bands());
	} else {
		while(meta.size() < tgtp.bands())
			meta.push_back("");
	}

	GridProps dstp(tgtp);
	dstp.setWritable(true);
	dstp.setBandMetadata(meta);

	geo::grid::Grid<float> dstg(dstfile, dstp);
	dstg.writeTo(tgtg);

	return 0;
}


