//============================================================================
// Name        : contrem.cpp
// Author      : Rob Skelly
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include "reader.hpp"
#include "writer.hpp"
#include "processor.hpp"

int main(int argc, char** argv) {

	int bufSize = 256;
	std::string datafile = argv[1];
	std::string bandfile = argv[2];
	std::string outfile = argv[3];

	GDALReader reader(datafile);
	reader.setBufSize(bufSize);

	GDALWriter writer(outfile, reader.cols(), reader.rows(), 1);

	Processor processor;

	std::vector<double> buf(bufSize * bufSize * reader.bands());
	std::vector<double> wavelengths = reader.getBands();
	std::vector<double> output(bufSize * bufSize * 1);

	int col, row, cols, rows;
	while(reader.next(buf, col, row, cols, rows)) {
		processor.process(buf, wavelengths, output, cols, rows, reader.bands(), bufSize);
		writer.write(output, col, row, cols, rows, bufSize);
	}

	return 0;
}
