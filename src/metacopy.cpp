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

	if(argc < 4) {
		std::cerr << "Usage: metacopy <meta source file> <target file> <property name>\n";
		return 1;
	}

	std::string srcfile = argv[1];
	std::string tgtfile = argv[2];
	std::string pname = argv[3];

	// Get the metadata from the source file.
	std::vector<std::string> meta;
	{
		geo::grid::Grid<float> srcg(srcfile);
		const GridProps& srcp = srcg.props();
		meta = srcp.bandMetadata();
	}

	// Initalize the target grid to add meta to.
	geo::grid::Grid<float> tgtg(tgtfile, true);
	tgtg.setMetadata(pname, meta);

	return 0;
}


