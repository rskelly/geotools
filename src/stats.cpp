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
inline void stats1(const std::vector<double>& v, double& n, double& min, double& max, double& mean) {
	min = SMAX;
	max = SMIN;
	n = 0;
	double sum = 0;
	for(double d : v) {
		if(d < min) min = d;
		if(d > max) max = d;
		sum += d;
		++n;
	}
	mean = sum / v.size();
}

#include <iomanip>
#include <fstream>

inline void stats2(const std::vector<double>& v, bool sample, double mean, double& variance, double& stddev, double& stderr, double& kurtosis, double& skewness, double& cov) {
	double sum2 = 0, sum3 = 0, sum4 = 0;
	for(double d : v) {
		sum2 += p2(d - mean);
		sum3 += p3(d - mean);
		sum4 += p4(d - mean);
	}
	{
		std::ofstream tmp("data/stats.txt");
		tmp << std::setprecision(12);
		for(const double& vv : v)
			tmp << vv << ",";
	}
	int n = sample ? v.size() - 1 : v.size();
	variance = sum2 / n;
	stddev = std::sqrt(variance);
	stderr = stddev / std::sqrt(sample ? n + 1 : n);
	skewness = (sum3 / n) / std::pow(std::sqrt(sum2), 3.0 / 2.0);
	kurtosis = (sum4 / n) / p2(variance);
	cov = stddev / mean;
}

void quantiles(const std::vector<double>& values, size_t quantiles, std::vector<double>& output) {
	size_t n = values.size();
	if(n <= quantiles) {
		for(size_t i = 0; i < values.size(); ++i)
			output.push_back(SNaN);
	} else {
		double step = (double) n / (quantiles - 1);
		output.resize(quantiles - 1);
		int i0, i1;
		for(size_t i = 1; i <= quantiles - 1; ++i) {
			i0 = (int) i * step;
			i1 = (int) std::ceil(i * step);
			output[i - 1] = (values[i1] + values[i0]) / 2.0;

		}
	}
}

inline void stats3(const std::vector<double>& v, double& median, double& mode, double* deciles, double& p25, double& p75, double& iqr2575) {

	std::vector<double> v0(v.begin(), v.end());
	std::sort(v0.begin(), v0.end());

	size_t n = v0.size();

	{
		// Compute the mode.
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
	}

	{
		// Compute the deciles.
		std::vector<double> dec;
		quantiles(v0, 10, dec);
		for(size_t i = 0; i < dec.size(); ++i)
			deciles[i] = dec[i];
	}

	{
		// Compute quartiles, iqr and median.
		std::vector<double> perc;
		int i0 = (int) (double) n / 99;
		int i1 = (int) std::ceil((double) n / 99);
		p25 = (v0[i0 * 25] + v0[i1 * 25]) / 2.0;
		p75 = (v0[i0 * 75] + v0[i1 * 75]) / 2.0;

		iqr2575 = p75 - p25;

		i0 = (int) (double) n / 2;
		i1 = (int) std::ceil((double) n / 2);
		median = (v0[i0] + v0[i1]) / 2.0;
	}

}


std::vector<std::string> __stats = {"n", "min", "max", "mean", "variance", "stddev", "stderr", "kurtosis", "skewness", "cov", "median", "mode", "p25", "p75", "iqr2575", "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9"};

std::vector<std::string> Stats::getStatNames() const {
	return __stats;
}

int Stats::computeStats(const std::vector<double>& values, std::vector<double>& output, bool sample) const {

	output.resize(__stats.size());

	stats1(values, output[0], output[1], output[2], output[3]); // n, min, max, mean
	stats2(values, sample, output[3], output[4], output[5], output[6], output[7], output[8], output[9]); // mean (in), variance, stddev, stderr, kurtosis, skewness, cov);
	stats3(values, output[9], output[10], (output.data() + 11), output[21], output[22], output[23]); // median, mode, deciles, p25, p75, iqr2575);

	return output.size();

}
