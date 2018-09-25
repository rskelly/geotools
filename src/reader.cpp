/*
 * reader.cpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <chrono>
#include <list>

#include <gdal_priv.h>

#include "reader.hpp"

using namespace hlrg;

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

void Reader::setBandMap(const std::map<int, int>& map) {
	m_bandMap = map;
	m_minIdx = 1;
	m_maxIdx = map.size();
	if(m_maxIdx > m_minIdx) {
		m_minWl = std::next(m_bandMap.begin(), m_minIdx - 1)->first;
		m_maxWl = std::next(m_bandMap.begin(), m_maxIdx - 1)->first;
	} else {
		m_minIdx = 0;
		m_minWl = 0;
		m_maxWl = 0;
	}
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

std::vector<double> Reader::getWavelengths() const {
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

std::vector<std::string> Reader::getBandNames() const {
	return m_bandNames;
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
			const char* m = band->GetMetadataItem(name.c_str());
			if(m) {
				// The wavelength is scaled so that exact matches can occur.
				int wl = (int) (atof(m) * WL_SCALE);
				if(wl > 0)
					bandMap[wl] = i;
			}
			m = band->GetDescription();
			if(m)
				m_bandNames.push_back(m);
		}
		if((int) bandMap.size() <= m_bands)
			std::runtime_error("The band map is incomplete -- wavelengths could not be read for all layers.");
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

	std::fill(buf.begin(), buf.end(), 0);

	col = m_col;
	row = m_row;
	cols = std::min(m_bufSize, m_cols - m_col);
	rows = std::min(m_bufSize, m_rows - m_row);

	m_col += m_bufSize;

	double* data = (double*) buf.data();
	for(int i = m_minIdx; i <= m_maxIdx; ++i) {
		//std::cerr << "band " << i << "\n";
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

		std::stringstream sbuf(buf);
		while(std::getline(sbuf, item, ' ')) {
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

	std::fill(buf.begin(), buf.end(), 0);

	col = m_col;
	row = m_row;
	cols = std::min(m_bufSize, m_cols - m_col);
	rows = std::min(m_bufSize, m_rows - m_row);

	m_col += m_bufSize;

	for(int i = m_minIdx; i <= m_maxIdx; ++i) {
		for(int r = row; r < std::min(row + m_bufSize, m_rows); ++r) {
			for(int c = col; c < std::min(col + m_bufSize, m_cols); ++c) {
				long px = ((long) c << 32) | r;
				if(m_pixels.find(px) != m_pixels.end()) {
					int idx = (i - m_minIdx) * m_bufSize * m_bufSize + r * m_bufSize + c;
					buf[idx] = m_pixels[px].values[i];
				}
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

		std::stringstream bufs(buf);
		while(std::getline(bufs, item, ',')) {
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

		fields.clear();
	}
}

const std::map<int, int>& BandMapReader::bandMap() const {
	return m_bandMap;
}




FrameIndexReader::FrameIndexReader(const std::string& filename) {
	std::ifstream in(filename, std::ios::in);
	std::string frame, time;
	size_t rpos;
	// Skip the header.
	if(in.good())
		std::getline(in, frame, '\n');
	// Read the data.
	std::list<std::pair<int, long> > items;
	while(in.good()) {
		std::getline(in, frame, ' '); // TODO: Configurable delimiter.
		std::getline(in, time, '\n');
		if(time.empty() || frame.empty())
			continue;
		if((rpos = time.find('\r')) != std::string::npos)
			time.replace(rpos, 1, 0, 'x');
		items.push_back(std::make_pair(std::stoi(frame), std::stol(time)));
	}

	{
		auto it = items.begin();
		std::advance(it, items.size() / 2);
		int mframe = it->first;
		int mtime = it->second;
		m_frames = new BinTree<long, int>(mtime, mframe);
		m_times = new BinTree<int, long>(mframe, mtime);
	}

	for(auto p : items) {
		m_frames->add(p.second, p.first);
		m_times->add(p.first, p.second);
	}
};

bool FrameIndexReader::getTime(int frame, long& time) const {
	long t;
	if(m_times->get(frame, t)) {
		time = t;
		return true;
	}
	return false;
}

bool FrameIndexReader::getNearestTime(int frame, int& actualFrame, long& time) const {
	long t;
	int f;
	if(m_times->findNearest(frame, f, t)) {
		actualFrame = f;
		time = t;
		return true;
	}
	return false;
}

bool FrameIndexReader::getFrame(long time, int& frame) const {
	int f;
	if(m_frames->get(time, f)) {
		frame = f;
		return true;
	}
	return false;
}

bool FrameIndexReader::getNearestFrame(long time, long& actualTime, int& frame) const {
	long t;
	int f;
	if(m_frames->findNearest(time, t, f)) {
		actualTime = t;
		frame = f;
		return true;
	}
	return false;
}

FrameIndexReader::~FrameIndexReader() {
	delete m_frames;
	delete m_times;
}

long _getUTCMilSec(const std::string& input, const std::string& fmt) {
	std::tm t = {};
	strptime(input.c_str(), fmt.c_str(), &t);
	double a = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::from_time_t(std::mktime(&t)).time_since_epoch()).count();
	std::string sfrac = input.substr(input.find('.'), 6);
	double b = std::stod(sfrac) * 1000;
	return a + b;
}

IMUGPSRow::IMUGPSRow(std::istream& in, double msOffset) :
	index(0) {
	std::string buf;
	std::getline(in, buf, '\t');
	roll = std::stod(buf);
	std::getline(in, buf, '\t');
	pitch = std::stod(buf);
	std::getline(in, buf, '\t');
	yaw = std::stod(buf);
	std::getline(in, buf, '\t');
	lat = std::stod(buf);
	std::getline(in, buf, '\t');
	lon = std::stod(buf);
	std::getline(in, buf, '\t');
	alt = std::stod(buf);
	std::getline(in, buf, '\t');
	gpsTime = std::stol(buf);
	std::getline(in, buf, '\t');
	utcTime = _getUTCMilSec(buf, "%Y/%m/%d %H:%M:%S") + msOffset;
	std::getline(in, buf, '\t');
	status = std::stod(buf);
	std::getline(in, buf);
	heading = std::stod(buf);
}


IMUGPSReader::IMUGPSReader(const std::string& filename, double msOffset) :
	m_lastIndex(0) {
	m_in.open(filename, std::ios::in);
	// Skip the header.
	std::string buf;
	std::getline(m_in, buf, '\n');
	std::list<IMUGPSRow*> rows;
	while(std::getline(m_in, buf, '\n')) {
		try {
			rows.push_back(new IMUGPSRow(m_in, msOffset));
		} catch(...) {
			continue;
		}
	}
	{
		auto it = rows.begin();
		std::advance(it, rows.size() / 2);
		IMUGPSRow* tmp = *it;
		m_gpsTimesT = new BinTree<long, IMUGPSRow*>(tmp->gpsTime, tmp);
		m_utcTimesT = new BinTree<long, IMUGPSRow*>(tmp->utcTime, tmp);
	}
	size_t i = 0;
	long mint = std::numeric_limits<long>::max(), maxt = std::numeric_limits<long>::lowest();
	m_rows.resize(rows.size());
	for(IMUGPSRow* row : rows) {
		row->index = i++;
		m_rows[row->index] = row;
		m_gpsTimesT->add(row->gpsTime, row);
		m_utcTimesT->add(row->utcTime, row);
		if(row->utcTime < mint) mint = row->utcTime;
		if(row->utcTime > maxt) maxt = row->utcTime;
	}
	std::cerr << "IMUGPS min date: " << mint << ", max date: " << maxt << "\n";
}

bool IMUGPSReader::getUTCTime(long gpsTime, long& utcTime) {
	long gpsTime0;
	IMUGPSRow* row;
	if(m_gpsTimesT->findNearest(gpsTime, gpsTime0, row)) {
		if(gpsTime == gpsTime0) {
			utcTime = row->utcTime;
			return true;
		} else if(gpsTime > gpsTime0) {
			if(row->index >= m_rows.size()) {
				// Past the last row, return the last value.
				utcTime = row->utcTime;
			} else {
				IMUGPSRow* row0 = m_rows[row->index];
				IMUGPSRow* row1 = m_rows[row->index + 1];
				// Interpolate the utc time, linearly.
				utcTime = row0->utcTime + (row1->utcTime - row0->utcTime) * gpsTime / (row1->gpsTime - row0->gpsTime);
			}
			return true;
		} else if(gpsTime < gpsTime0) {
			if(row->index == 0) {
				// At or before the first row, return the first value.
				utcTime = row->utcTime;
			} else {
				IMUGPSRow* row0 = m_rows[row->index - 1];
				IMUGPSRow* row1 = m_rows[row->index];
				// Interpolate the utc time, linearly.
				utcTime = row0->utcTime + (row1->utcTime - row0->utcTime) * gpsTime / (row1->gpsTime - row0->gpsTime);
			}
			return true;
		}
	}
	return false;
}

bool IMUGPSReader::getGPSTime(long utcTime, long& gpsTime) {
	long utcTime0;
	IMUGPSRow* row;
	if(m_utcTimesT->findNearest(utcTime, utcTime0, row)) {
		if(utcTime == utcTime0) {
			gpsTime = row->gpsTime;
			return true;
		} else if(utcTime > utcTime0) {
			if(row->index >= m_rows.size() - 1) {
				gpsTime = row->gpsTime;
			} else {
				IMUGPSRow* row0 = m_rows[row->index];
				IMUGPSRow* row1 = m_rows[row->index + 1];
				gpsTime = row0->gpsTime + (row1->gpsTime - row0->gpsTime) * (utcTime - row0->utcTime) / (row1->utcTime - row0->utcTime);
			}
			return true;
		} else if(utcTime < utcTime0) {
			if(row->index == 0) {
				gpsTime = row->gpsTime;
			} else {
				IMUGPSRow* row0 = m_rows[row->index - 1];
				IMUGPSRow* row1 = m_rows[row->index];
				gpsTime = row0->gpsTime + (row1->gpsTime - row0->gpsTime) * (utcTime - row0->utcTime) / (row1->utcTime - row0->utcTime);
			}
			return true;
		}
	}
	return false;
}

IMUGPSReader::~IMUGPSReader() {
	delete m_gpsTimesT;
	delete m_utcTimesT;
}

bool FlameRow::read(std::istream& in, double msOffset) {
	std::string buf;
	if(!std::getline(in, buf, ','))
		return false;
	dateTime = _getUTCMilSec(buf, "%Y-%m-%d %H:%M:%S") + msOffset;
	if(!std::getline(in, buf, ','))
		return false;
	utcTime = std::stol(buf) + msOffset;
	if(!std::getline(in, buf, '\n'))
		return false;
	size_t i = 0;
	std::stringstream cols(buf);
	while(std::getline(cols, buf, ','))
		bands[i++] = std::stod(buf);
	return in.good();
}

FlameReader::FlameReader(const std::string& filename, double msOffset) :
	m_msOffset(msOffset) {
	m_in.open(filename, std::ios::in);
	std::string buf;
	std::getline(m_in, buf, ',');
	std::getline(m_in, buf, ',');
	std::getline(m_in, buf, '\n');
	std::stringstream cols(buf);
	while(std::getline(cols, buf, ','))
		wavelengths.push_back(std::stod(buf));
}

bool FlameReader::next(FlameRow& row) {
	if(row.wavelengths.empty()) {
		row.wavelengths.assign(wavelengths.begin(), wavelengths.end());
		row.bands.resize(row.wavelengths.size());
	}
	return row.read(m_in, m_msOffset);
}

