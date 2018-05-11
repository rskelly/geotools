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
	void process(Reader& reader, const std::string& outfile, int bufSize, int threads);
};

#endif /* PROCESSOR_HPP_ */
