/*
 * stats.cpp
 *
 *  Created on: May 11", " 2018
 *      Author: rob
 */

#include <limits>
#include <unordered_map>
#include <algorithm>
#include <cmath>

#include "stats.hpp"

#define SMIN std::numeric_limits<double>::lowest()
#define SMAX std::numeric_limits<double>::max()
#define SNaN std::numeric_limits<double>::quiet_NaN()

#define p2(x) ((x)*(x))
#define p3(x) ((x)*(x)*(x))
#define p4(x) ((x)*(x)*(x)*(x))

using namespace hlrg;


// TODO: Should I ignore zeroes?
inline void stats1(const std::vector<double>& v, Stats& stats) {
	stats.min = SMAX;
	stats.max = SMIN;
	stats.n = 0;
	double sum = 0;
	for(double d : v) {
		if(d < stats.min) stats.min = d;
		if(d > stats.max) stats.max = d;
		sum += d;
	}
	stats.n = v.size();
	stats.mean = sum / stats.n;
}

#include <iomanip>
#include <fstream>

inline void stats2(const std::vector<double>& v, Stats& stats, bool sample) {
	double sum2 = 0, sum3 = 0, sum4 = 0;
	for(double d : v) {
		sum2 += p2(d - stats.mean);
		sum3 += p3(d - stats.mean);
		sum4 += p4(d - stats.mean);
	}
	int n = sample ? v.size() - 1 : v.size();
	stats.variance = sum2 / n;
	stats.stddev = std::sqrt(stats.variance);
	stats.stderr = stats.stddev / std::sqrt(sample ? n + 1 : n);
	stats.skewness = (sum3 / n) / std::pow(std::sqrt(sum2), 3.0 / 2.0);
	stats.kurtosis = (sum4 / n) / p2(stats.variance);
	stats.cov = stats.stddev / stats.mean;
}

void quantiles(const std::vector<double>& values, size_t quantiles, std::vector<double>& output) {
	size_t n = values.size();
	if(n <= quantiles) {
		for(size_t i = 1; i <  quantiles; ++i)
			output.push_back(SNaN);
	} else {
		for(size_t i = 1; i < quantiles; ++i) {
			size_t idx = (size_t) std::ceil((double) i / quantiles * n);
			output.push_back(values[idx]);
		}
	}
}

inline void stats3(const std::vector<double>& v, Stats& stats) {

	std::vector<double> v0(v.begin(), v.end());
	std::sort(v0.begin(), v0.end());

	{
		// Compute the mode.
		std::unordered_map<double, int> map;
		for(double d : v)
			map[d]++;
		int c = 1;
		stats.mode = SNaN;
		for(const auto& p : map) {
			if(p.second > c) {
				c = p.second;
				stats.mode = p.first;
			}
		}
	}

	{
		// Compute the deciles.
		std::vector<double> dec;
		quantiles(v0, 10, dec);
		stats.deciles.resize(9);
		for(size_t i = 0; i < dec.size(); ++i)
			stats.deciles[i] = dec[i];
	}

	stats.p25 = stats.p75 = stats.iqr = SNaN;
	if(v0.size() > 100) {
		// Compute quartiles, iqr and median.
		std::vector<double> perc;
		quantiles(v0, 100, perc);
		stats.p25 = perc[24];
		stats.p75 = perc[74];
		stats.iqr = stats.p75 - stats.p25;
	}

	stats.median = SNaN;
	if(!v0.empty()) {
		if(v0.size() % 2 == 0) {
			size_t n = v0.size() / 2;
			stats.median = (v0[n] + v0[n - 1]) / 2.0;
		} else {
			stats.median = v0[v0.size() / 2];
		}
	}
}


std::vector<std::string> __stats = {"n", "min", "max", "mean", "variance", "stddev", "stderr", "kurtosis", "skewness", "cov", "median", "mode", "p25", "p75", "iqr", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9"};

std::vector<std::string> Stats::getStatNames() {
	return __stats;
}

Stats Stats::computeStats(const std::vector<double>& values, bool sample) {

	Stats stats;

	stats1(values, stats); // n, min, max, mean
	stats2(values, stats, sample); // mean (in), variance, stddev, stderr, kurtosis, skewness, cov);
	stats3(values, stats); // median, mode, deciles, p25, p75, iqr2575);

	return stats;

}

std::vector<double> Stats::getStats() {
	std::vector<double> stats = {(double) n, min, max, mean, variance, stddev, stderr, kurtosis, skewness, cov, median, mode, p25, p75, iqr};
	for(double d : deciles)
		stats.push_back(d);
	return stats;
}

