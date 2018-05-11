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

// TODO: Should I ignore zeroes?
inline void stats1(const std::vector<double>& v, double& min, double& max, double& mean) {
	min = SMAX;
	max = SMIN;
	double sum = 0;
	for(double d : v) {
		if(d < min) min = d;
		if(d > max) max = d;
		sum += d;
	}
	mean = sum / v.size();
}

inline void stats2(const std::vector<double>& v, bool sample, double mean, double& variance, double& stddev, double& stderr, double& kurtosis, double& skewness, double& cov) {
	double sum2 = 0, sum3 = 0, sum4 = 0;
	for(double d : v) {
		sum2 += p2(d - mean);
		sum3 += p3(d - mean);
		sum4 += p4(d - mean);
	}
	int n = sample ? v.size() - 1 : v.size();
	variance = sum2 / n;
	stddev = std::sqrt(variance);
	stderr = stddev / std::sqrt(n + 1);
	skewness = sum3 / p3(std::sqrt(sum2));
	kurtosis = sum4 / p2(sum2);
	cov = stddev / mean;
}

inline void stats3(const std::vector<double>& v, double& median, double& mode, double* deciles, double& p25, double& p75, double& iqr2575) {

	std::vector<double> v0(v.begin(), v.end());
	std::sort(v0.begin(), v0.end());

	size_t n = v0.size();

	std::unordered_map<double, int> map;
	for(double d : v)
		map[d]++;

	int c = 1;
	mode = SNaN;
	for(const auto& p : map) {
		if(p.second > c) {
			c = p.second;
			mode = p.first;
		}
	}

	int step = n / 10;
	for(size_t i = 0; i < 10; ++i)
		deciles[i] = n >= 10 ? v0[(i + 1) * step] : SNaN;

	double idx = std::ceil(n / 100.0);
	p25 = v0[(int) idx * 25];
	p75 = v0[(int) idx * 75];
	iqr2575 = p75 - p25;

	idx = n / 2;
	median = n % 2 == 1 ? v[idx] : (v[idx - 1] + v[idx]) / 2.0;
}


std::vector<std::string> __stats = {"min", "max", "mean", "variance", "stddev", "stderr", "kurtosis", "skewness", "cov", "median", "mode", "p25", "p75", "iqr2575", "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9"};

std::vector<std::string> Stats::getStatNames() const {
	return __stats;
}

int Stats::computeStats(const std::vector<double>& values, std::vector<double>& output, bool sample) const {

	output.resize(__stats.size());

	stats1(values, output[0], output[1], output[2]); // min, max, mean
	stats2(values, sample, output[2], output[3], output[4], output[5], output[6], output[7], output[8]); // mean (in), variance, stddev, stderr, kurtosis, skewness, cov);
	stats3(values, output[8], output[9], (output.data() + 10), output[20], output[21], output[22]); // median, mode, deciles, p25, p75, iqr2575);

	return output.size();

}
