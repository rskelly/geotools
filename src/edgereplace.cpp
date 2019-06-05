/*
 * edgereplace.cpp
 *
 * This program removes the edges of a raw raster file up to a given width.
 *
 *  Created on: Aug 3, 2018
 *      Author: rob
 */

#include <string>
#include <fstream>
#include <vector>
#include <iostream>

inline void replaceFloat(float* buf, int imgWidth, int edgeWidth, float replace) {
	for(int i = 0; i < edgeWidth; ++i) {
		buf[i] = replace;
		buf[imgWidth - edgeWidth + i] = replace;
	}
}

inline void replaceUInt16(unsigned short* buf, int imgWidth, int edgeWidth, unsigned short replace) {
	for(int i = 0; i < edgeWidth; ++i) {
		buf[i] = replace;
		buf[imgWidth - edgeWidth + i] = replace;
	}
}

inline void trimFloat(float* buf, int imgWidth, int edgeWidth) {
	for(int i = edgeWidth; i < imgWidth - edgeWidth; ++i)
		buf[i - edgeWidth] = buf[i];
}

inline void trimUInt16(unsigned short* buf, int imgWidth, int edgeWidth) {
	for(int i = edgeWidth; i < imgWidth - edgeWidth; ++i)
		buf[i - edgeWidth] = buf[i];
}

int main(int argc, char** argv) {

	std::string infile = argv[1];
	std::string outfile = argv[2];
	int imgWidth = atoi(argv[3]);
	int edgeWidth = atoi(argv[4]);
	std::string replace = argv[5];
	int dataType = atoi(argv[6]);

	int dataSize = 2;

	switch(dataType) {
	case 4:
		dataSize = 4;
		break;
	case 12:
		dataSize = 2;
		break;
	}

	std::ifstream input(infile, std::ios::binary);
	std::ofstream output(outfile, std::ios::binary);

	std::vector<char> buf(dataSize * imgWidth);

	while(input.read(buf.data(), dataSize * imgWidth)) {
		switch(dataType) {
		case 4:
			//replaceFloat((float*) buf.data(), imgWidth, edgeWidth, (float) atof(replace.c_str()));
			trimFloat((float*) buf.data(), imgWidth, edgeWidth);
			break;
		case 12:
			//replaceUInt16((unsigned short*) buf.data(), imgWidth, edgeWidth, (unsigned short) atoi(replace.c_str()));
			trimUInt16((unsigned short*) buf.data(), imgWidth, edgeWidth);
			break;
		}
		output.write(buf.data(), dataSize * (imgWidth - edgeWidth * 2));
	}

	std::cerr << "Done.\n";
}

