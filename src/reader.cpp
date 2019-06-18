/*
 * reader.cpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#include <sys/mman.h>
#include <fcntl.h>
#include <sys/resource.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <chrono>
#include <list>
#include <iomanip>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include "reader.hpp"
#include "util.hpp"

using namespace hlrg::reader;
using namespace hlrg::ds;
using namespace hlrg::util;

namespace {

	long getUTCMilSec(const std::string& input, const std::string& fmt) {
		// hack: https://stackoverflow.com/questions/14504870/convert-stdchronotime-point-to-unix-timestamp#14505248
		std::tm t = {};
		std::stringstream ss(input);
		ss >> std::get_time(&t, fmt.c_str());
		time_t t1 = std::mktime(&t);
		t = *std::gmtime(&t1);
		time_t t2 = std::mktime(&t);
		double a = (t1 - (t2 - t1) + std::stod(input.substr(input.find('.'), 6))) * 1000;
		return a;
	}

}

Reader::Reader() :
	m_cols(0), m_rows(0), m_bands(0),
	m_col(0), m_row(0),
	m_bufSize(256),
	m_minWl(0), m_maxWl(0),
	m_minIdx(0), m_maxIdx(0) {
}

FileType Reader::fileType() const {
	return getFileType(m_filename);
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

std::map<int, int> Reader::getBandMap() {
	return m_bandMap;
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
	std::vector<std::string> names;
	if(m_maxIdx > 0 && m_maxIdx > m_minIdx) {
		for(auto p = std::next(m_bandNames.begin(), m_minIdx - 1); p != std::next(m_bandNames.begin(), m_maxIdx); ++p)
			names.push_back(*p);
		return names;
	}
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


PointSetReader::PointSetReader(const std::string& filename, const std::string& layer) :
	m_tree(nullptr) {

	GDALAllRegister();

	GDALDataset* ds = (GDALDataset*) GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR|GDAL_OF_READONLY, 0, 0, 0);
	if(!ds)
		throw std::runtime_error("Failed to open dataset " + filename);
	OGRLayer* lyr = ds->GetLayerByName(layer.c_str());
	if(!lyr)
		throw std::runtime_error("No layer named " + layer + " on dataset.");

	m_tree = new geo::ds::KDTree<hlrg::reader::Point>(2);

	OGRwkbGeometryType type = lyr->GetGeomType();
	if(type != wkbPoint && type != wkbMultiPoint)
		throw std::runtime_error("Illegal geometry type: " + std::to_string(type));

	OGRFeature* feat;
	OGRPoint* pt;
	OGRMultiPoint* mpt;
	while((feat = lyr->GetNextFeature())) {
		switch(type) {
		case wkbPoint:
			pt = (OGRPoint*) feat->GetGeometryRef();
			m_tree->add(new hlrg::reader::Point(pt->getX(), pt->getY(), 0));
			break;
		case wkbMultiPoint:
			mpt = (OGRMultiPoint*) feat->GetGeometryRef();
			for(int i = 0; i < mpt->getNumGeometries(); ++i) {
				pt = (OGRPoint*) mpt->getGeometryRef(i);
				m_tree->add(new hlrg::reader::Point(pt->getX(), pt->getY(), 0));
			}
			break;
		default:
			break;
		}
		OGRFeature::DestroyFeature(feat);
	}

	GDALClose(ds);

	m_tree->build();
}

std::vector<std::string> PointSetReader::getLayerNames(const std::string& filename) {

	GDALAllRegister();

	std::vector<std::string> names;

	void* h = GDALOpenEx(filename.c_str(), GDAL_OF_VECTOR, 0, 0, 0);
	GDALDataset* ds = (GDALDataset*) h;
	if(!ds)
		throw std::runtime_error("Failed to open dataset " + filename);
	for(int i = 0; i < ds->GetLayerCount(); ++i) {
		OGRLayer* lyr = ds->GetLayer(i);
		names.emplace_back(lyr->GetName());
	}

	return names;
}

int PointSetReader::search(double x, double y, double radius, std::vector<hlrg::reader::Point*>& pts) {

	std::vector<double> dist;
	return m_tree->radSearch(hlrg::reader::Point(x, y), radius, 0, std::back_inserter(pts), std::back_inserter(dist));
}

int PointSetReader::samplesNear(double x, double y, double radius) {

	std::vector<double> dist;
	std::vector<hlrg::reader::Point*> pts;
	return m_tree->radSearch(hlrg::reader::Point(x, y), radius, 0, std::back_inserter(pts), std::back_inserter(dist));
}

PointSetReader::~PointSetReader() {
	if(m_tree)
		delete m_tree;
}



GDALReader::GDALReader(const std::string& filename) : Reader(),
		m_ds(nullptr),
		m_mappedSize(0),
		m_mapped(nullptr) {

	GDALAllRegister();

	if(!(m_ds = (GDALDataset*) GDALOpen(filename.c_str(), GA_ReadOnly)))
		throw std::invalid_argument("Failed to open dataset.");

	m_minIdx = m_maxIdx = 1;
	m_bands = m_ds->GetRasterCount();
	m_cols = m_ds->GetRasterXSize();
	m_rows = m_ds->GetRasterYSize();
	m_ds->GetGeoTransform(m_trans);

	loadBandMap();
}

double GDALReader::toX(int col) const {
	return m_trans[0] + col * m_trans[1] + m_trans[1] * 0.5;
}

double GDALReader::toY(int row) const {
	return m_trans[3] + row + m_trans[5] + m_trans[5] * 0.5;
}

double GDALReader::toCol(double x) const {
	return (x - m_trans[0]) / m_trans[1];
}

double GDALReader::toRow(double y) const {
	return (y - m_trans[3]) / m_trans[5];
}

double GDALReader::resX() const {
	return m_trans[1];
}

double GDALReader::resY() const {
	return m_trans[5];
}

void GDALReader::remap() {
	remap(1, m_bands);
}

void GDALReader::remap(double minWl, double maxWl) {
	int iminWl = (int) std::round(minWl * WL_SCALE);
	int imaxWl = (int) std::round(maxWl * WL_SCALE);
	remap(m_bandMap[iminWl], m_bandMap[imaxWl]);
}

template <class T>
void doRemap(GDALDataset* ds, double* mapped, int minBand, int maxBand, int cols, int rows) {

	int bcols, brows, acols, arows;
	GDALRasterBand* firstBand = ds->GetRasterBand(minBand);
	int mappedBands = maxBand - minBand + 1;

	firstBand->GetBlockSize(&bcols, &brows);

	std::vector<T> buf(mappedBands * bcols * brows * sizeof(T));
	std::vector<double> row(mappedBands);

	for(int br = 0; br < rows / brows; ++br) {
		for(int bc = 0; bc < cols / bcols; ++bc) {
			std::cout << "Remapping block " << (br * (cols / bcols) + bc) << " of " << bcols * brows << "\n";

			// Get a "stack" of blocks representing the band data within a region of pixels.
			// This is BSQ oriented.
			for(int i = minBand, b = 0; i <= maxBand; ++i, ++b) {
				GDALRasterBand* band = ds->GetRasterBand(i);
				if(CE_None != band->ReadBlock(bc, br, buf.data() + (b * bcols * brows)))
					continue;
			}

			firstBand->GetActualBlockSize(bc, br, &acols, &arows);

			for(int r = 0; r < arows; ++r) {
				for(int c = 0; c < acols; ++c) {

					// Copy the band values into the row buffer.
					for(int b = 0; b < mappedBands; ++b)
						row[b] = buf[(b * bcols * brows) + r * bcols + c];

					// Write the row into the mapped file as a sequence of band values for the pixel.
					size_t idx = (br * brows + r) * cols * mappedBands + (bc * bcols + c) * mappedBands;
					std::memcpy(mapped + idx, row.data(), row.size() * sizeof(double));
				}
			}
		}
	}

}

void GDALReader::remap(int minBand, int maxBand) {
	m_mappedMinBand = minBand;
	m_mappedBands = (maxBand - minBand) + 1;
	m_mappedSize = (size_t) m_cols * m_rows * m_mappedBands * sizeof(double);
	m_mappedFile.reset(new TmpFile(m_mappedSize));
	m_mapped = (double*) mmap(0, m_mappedSize, PROT_READ|PROT_WRITE, MAP_SHARED, m_mappedFile->fd, 0);
	m_mappedFile->close();

	if((long) m_mapped == -1)
		throw std::runtime_error(std::string("Failed to remap: ") + strerror(errno) + " " + std::to_string(errno));

	GDALRasterBand* firstBand = m_ds->GetRasterBand(minBand);
	GDALDataType type = firstBand->GetRasterDataType();

	switch(type) {
	case GDT_Float32:
		doRemap<float>(m_ds, m_mapped, minBand, maxBand, cols(), rows());
		break;
	case GDT_Float64:
		doRemap<double>(m_ds, m_mapped, minBand, maxBand, cols(), rows());
		break;
	default:
		throw std::runtime_error("remap only implemented for float rasters.");
	}
}

double GDALReader::mapped(int col, int row, double wl) {
	return mapped(col, row, m_bandMap[(int) (wl * WL_SCALE)]);
}

double GDALReader::mapped(int col, int row, int band) {
	size_t idx = row * m_cols * m_mappedBands + col * m_mappedBands + (band - m_mappedMinBand);
	if(idx >= m_mappedSize)
		return std::nan("");
	return m_mapped[idx];
}

bool GDALReader::mapped(int col, int row, std::vector<double>& values) {
	size_t idx = row * m_cols * m_mappedBands + col * m_mappedBands;
	if(idx + m_mappedBands > m_mappedSize)
		return false;
	values.resize(m_mappedBands);
	std::memcpy(values.data(), m_mapped + idx, m_mappedBands * sizeof(double));
	return true;
}


void GDALReader::loadBandMap() {
	std::map<int, int> bandMap;
	std::vector<std::string> names = {"wavelength", "WAVELENGTH", "Description"};
	for(int i = 1; i <= m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		const char* m = nullptr;
		for(const std::string& name : names) {
			if((m = band->GetMetadataItem(name.c_str())) != nullptr)
				break;
		}
		if(m) {
			// The wavelength is scaled so that exact matches can occur.
			int wl = (int) (atof(m) * WL_SCALE);
			if(wl > 0)
				bandMap[wl] = i;
		}
		m = band->GetDescription();
		if(m) {
			m_bandNames.push_back(m);
			double wl = atof(m);
			if(wl > 0 && !std::isnan(wl) && bandMap.find(wl) == bandMap.end()) {
				// The wavelength is scaled so that exact matches can occur.
				int iwl = (int) (wl * WL_SCALE);
				bandMap[iwl] = i;
			}
		}
	}
	if((int) bandMap.size() < m_bands)
		std::runtime_error("The band map is incomplete -- wavelengths could not be read for all layers.");

	setBandMap(bandMap);
}

int GDALReader::toCol(double x) {
	double trans[6];
	m_ds->GetGeoTransform(trans);
	return (int) (x - trans[0]) / trans[1];
}

int GDALReader::toRow(double y) {
	double trans[6];
	m_ds->GetGeoTransform(trans);
	return (int) (y - trans[3]) / trans[5];
}

int GDALReader::getInt(double x, double y) {
	return GDALReader::getInt(toCol(x), toRow(y));
}

int GDALReader::getInt(int col, int row) {
	return (int) GDALReader::getFloat(col, row);
}

float GDALReader::getFloat(double x, double y) {
	return GDALReader::getInt(toCol(x), toRow(y));
}

float GDALReader::getFloat(int col, int row) {
	if(col < 0 || col >= m_ds->GetRasterXSize() || row < 0 || row >= m_ds->GetRasterYSize())
		return 0;
	float buf[1];
	GDALRasterBand* band = m_ds->GetRasterBand(1);
	if(CE_None != band->RasterIO(GF_Read, col, row, 1, 1, buf, 1, 1, GDT_Float32, 0, 0, 0))
		return 0;
	return buf[0];
}

bool GDALReader::next(std::string& id, std::vector<double>& buf, int& cols, int& col, int& row) {

	id = "";
	cols = m_cols;
	col = m_col;
	row = m_row;

	if(m_row >= m_rows)
		return false;

	if(m_mapped) {

		if(!mapped(m_col, m_row, buf))
			return false;

		if(++m_col >= m_cols) {
			m_col = 0;
			++m_row;
		}

	} else {

		int numBands = m_maxIdx - m_minIdx + 1;

		buf.resize(m_cols * numBands);
		std::fill(buf.begin(), buf.end(), 0);

		double* data = (double*) buf.data();
		for(int i = m_minIdx; i <= m_maxIdx; ++i) {
			//std::cerr << "band " << i << "\n";
			GDALRasterBand* band = m_ds->GetRasterBand(i);
			if(CE_None != band->RasterIO(GF_Read, 0, m_row, m_cols, 1, (void*) (data + i - m_minIdx), m_cols, 1, GDT_Float64, numBands * sizeof(double), 0, 0))
				return false;
		}

		++m_row;
	}

	return true;
}

bool GDALReader::next(std::vector<double>& buf, int band, int& cols, int& col, int& row) {

	if(m_row >= m_rows)
		return false;

	cols = m_cols;
	col = m_col;
	row = m_row;

	if(m_mapped) {

		std::vector<double> _buf;
		if(!mapped(m_col, m_row, _buf))
			return false;
		buf.resize(1);
		buf[0] = _buf[0];

		if(++m_col >= m_cols) {
			m_col = 0;
			++m_row;
		}

	} else {

		buf.resize(m_cols);
		std::fill(buf.begin(), buf.end(), 0);

		GDALRasterBand* bd = m_ds->GetRasterBand(band);
		if(CE_None != bd->RasterIO(GF_Read, 0, m_row, m_cols, 1, buf.data(), m_cols, 1, GDT_Float64, 0, 0, 0))
			return false;

		++m_row;
	}

	return true;
}

GDALReader::~GDALReader() {
	GDALClose(m_ds);
	if(m_mapped)
		munmap(m_mapped, m_mappedSize);
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
		std::getline(in, frame, '\t'); // TODO: Configurable delimiter.
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
	utcTime = getUTCMilSec(buf, "%Y/%m/%d %H:%M:%S") + msOffset;
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
		// Get the mid point of the list and create the bin trees using the mid values as root.
		auto it = rows.begin();
		std::advance(it, rows.size() / 2);
		IMUGPSRow* tmp = *it;
		m_gpsTimesT = new BinTree<long, IMUGPSRow*>(tmp->gpsTime, tmp);
		m_utcTimesT = new BinTree<long, IMUGPSRow*>(tmp->utcTime, tmp);
	}
	// Iterate over the rows, adding to the trees. The root values will just update.
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
	dateTime = getUTCMilSec(buf, "%Y-%m-%d %H:%M:%S") + msOffset;
	if(!std::getline(in, buf, ','))
		return false;
	utcTime = std::stol(buf);
	if(!std::getline(in, buf, '\n'))
		return false;
	size_t i = 0;
	std::stringstream cols(buf);
	while(std::getline(cols, buf, ','))
		bands[i++] = std::stod(buf);
	return in.good();
}

FlameReader::FlameReader(const std::string& filename, double msOffset) :
	m_msOffset(msOffset),
	m_filename(filename) {
	m_in.open(filename, std::ios::in);
	std::string buf;
	std::getline(m_in, buf, ',');
	std::getline(m_in, buf, ',');
	std::getline(m_in, buf, '\n');
	std::stringstream cols(buf);
	while(std::getline(cols, buf, ','))
		wavelengths.push_back(std::stod(buf));
}

int FlameReader::rows() {
	std::ifstream in(m_filename, std::ios::in);
	std::string buf;
	int count = 0;
	while(std::getline(in, buf, '\n'))
		++count;
	return count - 1; // -1 for the header row.
}

bool FlameReader::next(FlameRow& row) {
	if(row.wavelengths.empty()) {
		row.wavelengths.assign(wavelengths.begin(), wavelengths.end());
		row.bands.resize(row.wavelengths.size());
	}
	return row.read(m_in, m_msOffset);
}


CSVReader::CSVReader(const std::string& filename, bool transpose, int headerRows, int minWlCol, int maxWlCol, int idCol) :
	m_filename(filename),
	m_idx(0),
	m_transpose(transpose),
	m_minWlCol(minWlCol), m_maxWlCol(maxWlCol),
	m_headerRows(headerRows),
	m_idCol(idCol) {
	load();
}

void CSVReader::load() {
	m_data.clear();
	std::ifstream input(m_filename);
	m_cols = 1;
	m_rows = 0;
	int maxCols = 0;
	std::string line;
	while(std::getline(input, line)) {
		if(line.empty())
			continue;
		std::vector<std::string> row;
		std::string field;
		std::stringstream ss(line);
		while(std::getline(ss, field, ','))
			row.push_back(field);
		m_data.push_back(row);
		++m_rows;
		maxCols = std::max(maxCols, (int) row.size());
	}

	for(std::vector<std::string>& row : m_data)
		row.resize(maxCols);

	if(m_transpose)
		transpose();

	loadBandMap();
	reset();
}

void CSVReader::loadBandMap() {
	std::map<int, int> map;
	size_t hIdx = m_headerRows > 0 ? m_headerRows - 1 : 0;
	if(m_data.size() > hIdx) {
		const std::vector<std::string>& header = m_data[hIdx];
		if(m_minWlCol >= 0 && m_maxWlCol < (int) header.size()) {
			for(int i = m_minWlCol; i <= m_maxWlCol; ++i) {
				double wl = atof(header[i].c_str());
				int idx = (int) (wl * WL_SCALE);
				map[idx] = i;
			}
		}
	}
	setBandMap(map);
}

void CSVReader::reset() {
	m_idx = m_headerRows;
}

void CSVReader::transpose() {
	std::vector<std::vector<std::string> > tmp(m_cols);
	for(int c = 0; c < m_cols; ++c)
		tmp[c].resize(m_rows);
	for(int r = 0; r < m_rows; ++r) {
		for(int c = 0; c < m_cols; ++c)
			tmp[c][r] = m_data[r][c];
	}
	m_data.swap(tmp);
	m_cols = m_rows;
	m_rows = m_data.size();
}

bool CSVReader::next(std::string& id, std::vector<double>& buf, int& cols, int& col, int& row) {
	if(m_idx >= m_rows)
		return false;

	cols = 1;
	col = 0;
	row = m_idx;

	std::vector<double> data(m_maxIdx - m_minIdx + 1);
	for(int i = m_minIdx; i <= m_maxIdx; ++i)
		data[i - m_minIdx] = atof(m_data[m_idx][i].c_str());
	buf.assign(data.begin(), data.end());

	id = m_data[m_idx][m_idCol];

	++m_idx;

	return true;
}


void CSVReader::guessFileProperties(const std::string& filename, bool& transpose, int& header, int& minCol, int& maxCol, int& idCol) {

	idCol = -1;
	header = 1;
	transpose = false;

	std::string line;
	std::string field;
	std::stringstream linestr;
	std::vector<bool> headerisfloat;
	std::vector<bool> isfloat;

	// Run through the first header/rows.
	{
		std::ifstream input(filename);
		int row = 0;
		while(std::getline(input, line)) {
			linestr.clear();
			linestr << line;
			if(row == header - 1) {
				while(std::getline(linestr, field, ','))
					headerisfloat.push_back(!field.empty() && atof(field.c_str()) != 0);
			} else if(row > header - 1) {
				while(std::getline(linestr, field, ','))
					isfloat.push_back(!field.empty() && atof(field.c_str()) != 0);
				break;
			}
			++row;
		}
	}

	// Find the first float col.
	minCol = 0;
	while(minCol < (int) headerisfloat.size() && !headerisfloat[minCol])
		++minCol;

	// Find he last float column.
	maxCol = (int) (headerisfloat.size() - 1);
	while(maxCol > minCol && !headerisfloat[maxCol])
		--maxCol;

	// Find the id column.
	for(size_t i = 0; i < isfloat.size(); ++i) {
		if(!isfloat[i]) {
			idCol = i;
			break;
		}
	}

	// Look for mismatches.
	if(isfloat.size() < headerisfloat.size())
		isfloat.resize(headerisfloat.size());
	for(int i = minCol; i <= maxCol; ++i) {
		if(headerisfloat[i] != isfloat[i]) {
			transpose = true;
			break;
		}
	}

	// If transpose not required, return.
	if(!transpose)
		return;

	{
		isfloat.clear();
		std::ifstream input(filename);
		int row = 0;
		while(std::getline(input, line)) {
			linestr.clear();
			linestr << line;
			std::getline(linestr, field, ',');
			isfloat.push_back(!field.empty() && atof(field.c_str()) != 0);
			++row;
		}
	}

	// Find the first float col.
	minCol = 0;
	while(minCol < (int) headerisfloat.size() && !headerisfloat[minCol])
		++minCol;

	// Find he last float column.
	maxCol = (int) (headerisfloat.size() - 1);
	while(maxCol > minCol && !headerisfloat[maxCol])
		--maxCol;
}
