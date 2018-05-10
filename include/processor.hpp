/*
 * processor.hpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#ifndef PROCESSOR_HPP_
#define PROCESSOR_HPP_

#include <vector>

class Processor {
public:
	void process(const std::vector<double>& buffer, const std::vector<double>& wavelengths,
			std::vector<double>& output, int cols, int rows, int bands, int bufSize);
};



#endif /* PROCESSOR_HPP_ */
