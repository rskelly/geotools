/*
 * processor.hpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#ifndef PROCESSOR_HPP_
#define PROCESSOR_HPP_

#include "reader.hpp"

class ProcessorConfig {
public:
	std::string outfile;
	std::string driver;
	std::string extension;
	int bufferSize;
	int threads;
	bool sampleStats;
	bool flagsRaster;
	double interpDistance;
};

class Processor {
public:
	void process(Reader* reader, const ProcessorConfig& config);
};

#endif /* PROCESSOR_HPP_ */
