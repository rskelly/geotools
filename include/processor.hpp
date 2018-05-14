/*
 * processor.hpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#ifndef PROCESSOR_HPP_
#define PROCESSOR_HPP_

#include "reader.hpp"

class Processor {
public:
	void process(Reader* reader, const std::string& outfile, const std::string& outDriver, int bufSize, int threads, bool sample);
};

#endif /* PROCESSOR_HPP_ */
