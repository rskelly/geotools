/*
 * reader.cpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#include <iostream>
#include <gdal_priv.h>

#include "reader.hpp"

GDALReader::GDALReader(const std::string& filename) :
	m_ds(nullptr),
	m_cols(0), m_rows(0), m_bands(0),
	m_col(0), m_row(0),
	m_bufSize(256),
	m_minWl(0), m_maxWl(0),
	m_minIdx(0), m_maxIdx(0) {

	GDALAllRegister();

	if(!(m_ds = (GDALDataset*) GDALOpen(filename.c_str(), GA_ReadOnly)))
		throw std::invalid_argument("Failed to open dataset.");

	m_bands = m_ds->GetRasterCount();
	m_cols = m_ds->GetRasterXSize();
	m_rows = m_ds->GetRasterYSize();

	{
		for(int i = 1; i <= m_bands; ++i) {
			GDALRasterBand* band = m_ds->GetRasterBand(i);
			char** meta = band->GetMetadata();
			if(*meta) {
				double wl = atof(*meta);
				if(wl > 0.0) {
					m_bandMap[wl] = i;
				} else {
					m_bandMap.clear();
					std::cerr << "Failed to read band map from metadata." << std::endl;
					break;
				}
			}
		}
	}
}

void GDALReader::setBufSize(int bufSize) {
	m_bufSize = bufSize;
}

bool GDALReader::next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows) {
	if(m_col >= m_cols) {
		m_row += m_bufSize;
		m_col = 0;
	}
	if(m_row >= m_rows)
		return false;
	col = m_col;
	row = m_row;
	cols = std::min(m_bufSize, m_cols - m_col);
	rows = std::min(m_bufSize, m_rows - m_row);
	double* data = (double*) buf.data();
	for(int i = 1; i <= m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Read, m_col, m_row, cols, rows,
				(void*) (data + (i - 1) * m_bufSize * m_bufSize),
				m_bufSize, m_bufSize, GDT_Float64, 0, 0, 0))
			return false;
	}
	return true;
}

void GDALReader::setBandMap(std::map<double, int>& map) {
	m_bandMap = map;
	m_minIdx = 0;
	m_maxIdx = map.size() - 1;
	m_minWl = m_bandMap.begin()->first;
	m_maxWl = std::next(m_bandMap.begin(), m_bandMap.size())->first;
}

void GDALReader::setBandRange(double min, double max) {
	for(const auto& pair : m_bandMap) {
		if(pair.first <= min)
			m_minIdx = std::max(0, pair.second - 1);
		if(pair.first >= max) {
			m_maxIdx = pair.second;
			break;
		}
	}
	m_minWl = min;
	m_maxWl = max;
}

std::vector<double> GDALReader::getBands() const {
	return {m_minWl, m_maxWl};
}

std::vector<int> GDALReader::getIndices() const {
	return {m_minIdx, m_maxIdx};
}

int GDALReader::cols() const {
	return m_cols;
}

int GDALReader::rows() const {
	return m_rows;
}

int GDALReader::bands() const {
	return m_bands;
}

GDALReader::~GDALReader() {
	GDALClose(m_ds);
}
