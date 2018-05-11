/*
 * stats.hpp
 *
 *  Created on: May 11, 2018
 *      Author: rob
 */

#ifndef INCLUDE_STATS_HPP_
#define INCLUDE_STATS_HPP_

#include <vector>
#include <string>

class Stats {
public:
	std::vector<std::string> getStatNames() const;
	/**
	 * Compute stats.
	 * If sample is true, use sample statistics. If false, use population.
	 */
	int computeStats(const std::vector<double>& values, std::vector<double>& output, bool sample = true) const;
};



#endif /* INCLUDE_STATS_HPP_ */
