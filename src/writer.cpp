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
#include <cstdlib>
#include <cstring>

#include "writer.hpp"
#include "stats.hpp"

bool __isnonzero(const double& v) {
	return v > 0;
}


int makedir(const std::string& filename) {
	std::string path = filename.substr(0, filename.find_last_of('/'));
	if(mkdir(path.c_str(), 0755) == 0)
		return 0;
	switch(errno) {
	case EEXIST: return 0;
	default: return errno;
	}
}

GDALWriter::GDALWriter(const std::string& filename, const std::string& driver, int cols, int rows, int bands,
		const std::string& fieldName, const std::vector<std::string>& bandNames, DataType type) :
	m_ds(nullptr),
	m_bands(0), m_cols(0), m_rows(0) {

	int err;
	if((err = makedir(filename)))
		throw std::runtime_error("Could not create directory for " + filename + "; " + std::to_string(err));

	GDALDataType gtype;
	switch(type) {
	case DataType::Byte: gtype = GDT_Byte; break;
	case DataType::Int32: gtype = GDT_Int32; break;
	case DataType::Float32: gtype = GDT_Float32; break;
	default:
		throw std::invalid_argument("Invalid data type.");
	}

	GDALAllRegister();
	GDALDriverManager* gm = GetGDALDriverManager();
	GDALDriver* drv = gm->GetDriverByName(driver.c_str());
	if(!drv)
		throw std::runtime_error("Driver not found: " + driver);

	m_ds = drv->Create(filename.c_str(), cols, rows, bands, gtype, nullptr);
	if(!m_ds)
		throw std::runtime_error("Failed to open " + filename + " for writing");
	m_bands = m_ds->GetRasterCount();
	m_cols = m_ds->GetRasterXSize();
	m_rows = m_ds->GetRasterYSize();

	if(!bandNames.empty()) {
		for(int i = 1; i <= std::min((int) bandNames.size(), m_bands); ++i) {
			if(CE_None != m_ds->GetRasterBand(i)->SetMetadataItem(fieldName.c_str(), bandNames[i - 1].c_str()))
				std::cerr << "Failed to set metadata item " << fieldName << "\n";
		}
	}
}

bool GDALWriter::write(const std::vector<double>& buf, int col, int row, int cols, int rows, int bufSize) {
	if(col < 0 || col >= m_cols || col + cols > m_cols
			|| row < 0 || row >= m_rows || row + rows > m_rows)
		return false;
	const double* data = buf.data();
	for(int i = 1; i <= m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Write, col, row, cols, rows,
				(void*) (data + (i - 1) * bufSize * bufSize),
				bufSize, bufSize, GDT_Float64, 0, 0, 0))
			return false;
	}
	return true;
}

bool GDALWriter::write(const std::vector<int>& buf, int col, int row, int cols, int rows, int bufSize) {
	if(col < 0 || col >= m_cols || col + cols > m_cols
			|| row < 0 || row >= m_rows || row + rows > m_rows)
		return false;
	const int* data = buf.data();
	for(int i = 1; i <= m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Write, col, row, cols, rows,
				(void*) (data + (i - 1) * bufSize * bufSize),
				bufSize, bufSize, GDT_Int32, 0, 0, 0))
			return false;
	}
	return true;
}


bool GDALWriter::writeStats(const std::string& filename, const std::vector<std::string>& names) {

	if(!names.empty() && (int) names.size() != m_bands)
		throw std::invalid_argument("Band names must be the same size as the number of bands, or empty.");

	Stats stats;

	std::vector<double> buf(m_cols * m_rows);
	std::vector<std::string> statNames = stats.getStatNames();
	std::vector<double> results(statNames.size());

	std::ofstream out(filename, std::ios::out);
	out << std::setprecision(12) << "name";
	for(const std::string& name : statNames)
		out << "," << name;
	out << "\n";

	m_ds->FlushCache();

	for(int i = 1; i <= m_bands; ++i) {

		// Read an entire band into the buffer.
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Read, 0, 0, m_cols, m_rows, (void*) buf.data(), m_cols, m_rows, GDT_Float64, 0, 0, 0))
			return false;

		// Filter the values list to eliminate zeroes.
		std::vector<double> values;
		std::copy_if(buf.begin(), buf.end(), std::back_inserter(values), __isnonzero);

		if(!values.empty()) {
			Stats s = Stats::computeStats(values);
			out << names[i - 1];
			for(double v : s.getStats())
				out << "," << v;
			out << "\n";
		}
	}
	return true;
}

GDALWriter::~GDALWriter() {
	GDALClose(m_ds);
}

