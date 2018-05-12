/*
 * reader.cpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#include <iostream>
#include <fstream>
#include <sstream>

#include <gdal_priv.h>

#include "reader.hpp"

#define WL_SCALE 100000

Reader::Reader() :
	m_cols(0), m_rows(0), m_bands(0),
	m_col(0), m_row(0),
	m_bufSize(256),
	m_minWl(0), m_maxWl(0),
	m_minIdx(0), m_maxIdx(0) {
}

void Reader::setBufSize(int bufSize) {
	m_bufSize = bufSize;
}

void Reader::setBandMap(std::map<int, int>& map) {
	m_bandMap = map;
	m_minIdx = 0;
	m_maxIdx = map.size() - 1;
	m_minWl = std::next(m_bandMap.begin(), m_minIdx)->first;
	m_maxWl = std::next(m_bandMap.begin(), m_maxIdx)->first;
}

void Reader::setBandMap(std::string& bandfile) {
	throw std::runtime_error("Not implemented.");
}

void Reader::setBandRange(double min, double max) {
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

std::vector<double> Reader::getBandRange() const {
	return {(double) m_minWl / WL_SCALE, (double) m_maxWl / WL_SCALE};
}

std::vector<double> Reader::getBands() const {
	std::vector<double> bands;
	if(m_maxIdx > 0 && m_maxIdx > m_minIdx) {
		for(auto p = std::next(m_bandMap.begin(), m_minIdx - 1); p != std::next(m_bandMap.begin(), m_maxIdx); ++p)
			bands.push_back((double) p->first / WL_SCALE);
	} else {
		for(const auto& p : m_bandMap)
			bands.push_back((double) p.first / WL_SCALE);
	}
	return bands;
}

std::vector<int> Reader::getIndices() const {
	return {m_minIdx, m_maxIdx};
}

int Reader::cols() const {
	return m_cols;
}

int Reader::rows() const {
	return m_rows;
}

int Reader::bands() const {
	return m_maxIdx - m_minIdx + 1;
}


GDALReader::GDALReader(const std::string& filename) : Reader(),
		m_ds(nullptr) {

	GDALAllRegister();

	if(!(m_ds = (GDALDataset*) GDALOpen(filename.c_str(), GA_ReadOnly)))
		throw std::invalid_argument("Failed to open dataset.");

	m_bands = m_ds->GetRasterCount();
	m_cols = m_ds->GetRasterXSize();
	m_rows = m_ds->GetRasterYSize();

	{
		std::map<int, int> bandMap;
		std::string name = "wavelength";
		for(int i = 1; i <= m_bands; ++i) {
			GDALRasterBand* band = m_ds->GetRasterBand(i);
			std::string meta(band->GetMetadataItem(name.c_str()));
			if(!meta.empty()) {
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

	m_col += m_bufSize;

	double* data = (double*) buf.data();
	for(int i = m_minIdx; i <= m_maxIdx; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Read, col, row, cols, rows,
				(void*) (data + (i - m_minIdx) * m_bufSize * m_bufSize),
				m_bufSize, m_bufSize, GDT_Float64, 0, 0, 0))
			return false;
	}
	return true;
}

GDALReader::~GDALReader() {
	GDALClose(m_ds);
}



ROIReader::ROIReader(const std::string& filename) : Reader() {

	// Attempt to read in the ROI file.
	std::ifstream input(filename, std::ios::in);
	std::string buf, item;
	std::vector<std::string> fields;
	while(std::getline(input, buf)) {
		if(buf[0] == ';') continue;

		while(std::getline(std::stringstream(item), buf, ' ')) {
			if(!item.empty())
				fields.push_back(item);
		}

		int col = atoi(fields[1].c_str());
		int row = atoi(fields[2].c_str());

		// Allocate and retrieve the pixel; set its coordinates.
		px& p = m_pixels[((long) col << 32) | row];
		p.c = col;
		p.r = row;

		// Read the band values for the pixel.
		for(size_t i = 3; i < fields.size(); ++i)
			p.values.push_back(atof(fields[i].c_str()));

		if(col > m_cols)
			m_cols = col;
		if(row > m_rows)
			m_rows = row;
		m_bands = p.values.size();

		fields.clear();
	}

}

bool ROIReader::next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows) {
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

	m_col += m_bufSize;

	for(int i = m_minIdx; i <= m_maxIdx; ++i) {
		for(int r = row; r < std::min(row + m_bufSize, m_rows); ++r) {
			for(int c = col; c < std::min(col + m_bufSize, m_cols); ++c) {
				int idx = (i - m_minIdx) * m_bufSize * m_bufSize + r * m_bufSize + c;
				long px = ((long) col << 32) | row;
				buf[idx] = m_pixels[px].values[i];
			}
		}
	}
	return true;
}

ROIReader::~ROIReader() {
}


BandMapReader::BandMapReader(const std::string& filename, int wlCol, int idxCol, bool hasHeader) {

	std::ifstream input(filename, std::ios::in);
	std::string buf, item;

	if(hasHeader) {
		// Remove header, but fail if not found.
		if(!std::getline(input, buf))
			throw std::runtime_error("Failed to read header from file.");
	}

	std::vector<std::string> fields;
	while(std::getline(input, buf)) {

		while(std::getline(std::stringstream(item), buf, ' ')) {
			if(!item.empty())
				fields.push_back(item);
		}

		if(wlCol < 0 || wlCol >= (int) fields.size())
			std::invalid_argument("The wavelength column is invalid.");
		if(idxCol < 0 || idxCol >= (int) fields.size())
			std::invalid_argument("The band index column is invalid.");

		int wl = (int) atof(fields[wlCol].c_str()) * WL_SCALE;
		int idx = atoi(fields[idxCol].c_str());

		m_bandMap[wl] = idx;
	}
}

const std::map<int, int>& BandMapReader::bandMap() const {
	return m_bandMap;
}
