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

namespace hlrg {

/**
 * Compute a range of statistics on a list of given values.
 */
class Stats {
public:

	int n;				///<! The number of values.
	double min;			///<! The minimum value.
	double max;			///<! The maximum value.
	double mean;		///<! The arithmetic mean of values.
	double median;		///<! The median value.
	double variance;	///<! The (population/sample) variance of values.
	double stddev;		///<! The (population/sample) standard deviation of values.
	double stderr;		///<! The (population/sample) standard error of values.
	double cov;			///<! The (population/sample) coefficient of variance of values.
	double kurtosis;	///<! The kurtosis of the distribution of values.
	double skewness;	///<! The skewness of the distribution of values.
	double mode;		///<! The mode of values.
	double p25;			///<! The 25th percentile of values.
	double p75;			///<! The 75th percentile of values.
	double iqr;			///<! The interquartile range.
	std::vector<double> deciles;	///<! A list of the deciles.

	/**
	 * Get the names of available statistics as a vector.
	 *
	 * @return The names of available statistics as a vector.
	 */
	static std::vector<std::string> getStatNames();

	/**
	 * Return the computed statistics as a vector.
	 *
	 * @return The computed statistics as a vector.
	 */
	std::vector<double> getStats();

	/**
	 * Compute the statistics.
	 * If sample is true, use sample statistics. If false, use population.
	 *
	 * @param values The values to compute statistics on.
	 * @param sample If true, compute sample statistics, otherwise population statistics.
	 * @return A Stats object containing the results.
	 */
	static Stats computeStats(const std::vector<double>& values, bool sample = true);

};

} // hlrg


#endif /* INCLUDE_STATS_HPP_ */
