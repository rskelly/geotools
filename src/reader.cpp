/*
 * reader.cpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#include <iostream>
#include <gdal_priv.h>

#include "reader.hpp"

#define WL_SCALE 100000

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
		std::map<int, int> bandMap;
		for(int i = 1; i <= m_bands; ++i) {
			GDALRasterBand* band = m_ds->GetRasterBand(i);
			std::string meta(*band->GetMetadata());
			if(!meta.empty() && meta[0] == 'w') {
				meta = meta.substr(meta.find('=') + 1);
				// The wavelength is scaled so that exact matches can occur.
				int wl = (int) (atof(meta.c_str()) * WL_SCALE);
				if(wl > 0) {
					bandMap[wl] = i;
				} else {
					bandMap.clear();
					std::cerr << "Failed to read band map from metadata." << std::endl;
					break;
				}
			}
		}
		setBandMap(bandMap);
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
	for(int i = m_minIdx; i <= m_maxIdx; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Read, m_col, m_row, cols, rows,
				(void*) (data + (i - m_minIdx) * m_bufSize * m_bufSize),
				m_bufSize, m_bufSize, GDT_Float64, 0, 0, 0))
			return false;
	}
	return true;
}

void GDALReader::setBandMap(std::map<int, int>& map) {
	m_bandMap = map;
	m_minIdx = 0;
	m_maxIdx = map.size() - 1;
	m_minWl = std::next(m_bandMap.begin(), m_minIdx)->first;
	m_maxWl = std::next(m_bandMap.begin(), m_maxIdx)->first;
}

void GDALReader::setBandMap(std::string& bandfile) {
	throw std::runtime_error("Not implemented.");
}

void GDALReader::setBandRange(double min, double max) {
	int mins = (int) (min * WL_SCALE);
	int maxs = (int) (max * WL_SCALE);
	for(auto it = m_bandMap.begin(); it != m_bandMap.end(); ++it) {
		if(it->first <= mins) {
			m_minWl = it->first;
			m_minIdx = it->second;
		}
		if(it->first >= maxs) {
			m_maxWl = it->first;
			m_maxIdx = it->second;
			break;
		}
	}
}

std::vector<double> GDALReader::getBandRange() const {
	return {(double) m_minWl / WL_SCALE, (double) m_maxWl / WL_SCALE};
}

std::vector<double> GDALReader::getBands() const {
	std::vector<double> bands;
	for(auto p = std::next(m_bandMap.begin(), m_minIdx); p != std::next(m_bandMap.begin(), m_maxIdx + 1); ++p)
		bands.push_back((double) p->first / WL_SCALE);
	return bands;
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
	return m_maxIdx - m_minIdx + 1;
}

GDALReader::~GDALReader() {
	GDALClose(m_ds);
}
