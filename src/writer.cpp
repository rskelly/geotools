/*
 * writer.cpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#include <gdal_priv.h>
#include <errno.h>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include "writer.hpp"

int makedir(const std::string& filename) {
	std::string path = filename.substr(0, filename.find_last_of('/'));
	if(mkdir(path.c_str(), 0755) == 0)
		return true;
	switch(errno) {
	case EEXIST: return 0;
	default: return errno;
	}
}

GDALWriter::GDALWriter(const std::string& filename, int cols, int rows, int bands) :
	m_ds(nullptr),
	m_bands(0), m_cols(0), m_rows(0) {

	int err;
	if((err = makedir(filename)))
		throw std::runtime_error("Could not create directory for " + filename + "; " + std::to_string(err));

	GDALAllRegister();
	GDALDriverManager* gm = GetGDALDriverManager();
	GDALDriver* drv = gm->GetDriverByName("GTiff");
	if(!drv)
		throw std::runtime_error("GeoTiff driver not found.");
	m_ds = drv->Create(filename.c_str(), cols, rows, bands, GDT_Float32, nullptr);
	m_bands = m_ds->GetRasterCount();
	m_cols = m_ds->GetRasterXSize();
	m_rows = m_ds->GetRasterYSize();
}

bool GDALWriter::write(std::vector<double>& buf, int col, int row, int cols, int rows, int bufSize) {
	if(col < 0 || col >= m_cols || col + cols > m_cols
			|| row < 0 || row >= m_rows || row + rows > m_rows)
		return false;
	double* data = buf.data();
	for(int i = 1; i <= m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Write, col, row, cols, rows,
				(void*) (data + (i - 1) * bufSize * bufSize),
				bufSize, bufSize, GDT_Float64, 0, 0, 0))
			return false;
	}
	return true;
}

#define SMIN std::numeric_limits<double>::lowest()
#define SMAX std::numeric_limits<double>::max()
#define SNaN std::numeric_limits<double>::quiet_NaN()
#define p2(x) (x*x)
#define p3(x) (x*x*x)
#define p4(x) (x*x*x*x)

inline void stats1(std::vector<double>& v, double& min, double& max, double& mean) {
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

inline void stats2(std::vector<double>& v, double mean, double& variance, double& stddev, double& kurtosis, double& skewness, double& cov) {
	double sum2 = 0, sum3 = 0, sum4 = 0;
	for(double d : v) {
		sum2 += p2(d - mean);
		sum3 += p3(d - mean);
		sum4 += p4(d - mean);
	}
	int n = v.size() - 1;
	variance = sum2 / n;
	skewness = sum3 / n * p3(mean);
	kurtosis = sum4 / n * p4(mean);
	stddev = std::sqrt(variance);
	cov = stddev / mean;
}

inline void stats3(std::vector<double>& v, double& median, double& mode) {
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
	int idx = v.size() / 2;
	median = v.size() % 2 == 1 ? v[idx] : (v[idx - 1] + v[idx]) / 2.0;
}

inline void stats4(std::vector<double>& v, double* deciles, double& p25, double& p75, double& iqr2575) {
	std::sort(v.begin(), v.end());
	int step = (int) std::ceil(v.size() / 10.0);
	for(size_t i = step, j = 0; i < v.size(); i += step, ++j)
		deciles[j] = v[i];
	double idx = std::ceil(v.size() / 100.0);
	p25 = v[(int) idx * 25];
	p75 = v[(int) idx * 75];
	iqr2575 = p75 - p25;
}

bool GDALWriter::writeStats(const std::string& filename) {
	//double min, max, mean, variance, stddev, kurtosis, skewness, cov, median, mode;
	//double p25, p75, iqr2575;
	std::vector<double> stats(23);
	std::vector<double> buf(m_cols * m_rows);

	std::ofstream out(filename, std::ios::out);
	out << std::setprecision(12);
	out << "min,max,mean,variance,stddev,kurtosis,skewness,cov,median,mode,p25,p75,iqr2575,d0,d1,d2,d3,d4,d5,d6,d7,d8,d9\n";

	m_ds->FlushCache();

	for(int i = 1; i <= m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Read, 0, 0, m_cols, m_rows, (void*) buf.data(), m_cols, m_rows, GDT_Float64, 0, 0, 0))
			return false;
		stats1(buf, stats[0], stats[1], stats[2]); // min, max, mean
		stats2(buf, stats[2], stats[3], stats[4], stats[5], stats[6], stats[7]); // mean (in), variance, stddev, kurtosis, skewness, cov);
		stats3(buf, stats[8], stats[9]); // median, mode);
		stats4(buf, (stats.data() + 10), stats[20], stats[21], stats[22]); // deciles, p25, p75, iqr2575);
		for(double s: stats)
			out << s << ",";
		out << "\n";
	}
	return true;
}

GDALWriter::~GDALWriter() {
	GDALClose(m_ds);
}

