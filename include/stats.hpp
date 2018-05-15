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

	int n;
	double min;
	double max;
	double mean;
	double median;
	double variance;
	double stddev;
	double stderr;
	double cov;
	double kurtosis;
	double skewness;
	double mode;
	double p25;
	double p75;
	double iqr;
	std::vector<double> deciles;

	static std::vector<std::string> getStatNames();

	std::vector<double> getStats();

	/**
	 * Compute stats.
	 * If sample is true, use sample statistics. If false, use population.
	 */
	static Stats computeStats(const std::vector<double>& values, bool sample = true);

};



#endif /* INCLUDE_STATS_HPP_ */
