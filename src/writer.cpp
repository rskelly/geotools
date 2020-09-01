/*
 * writer.cpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#include <errno.h>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <sstream>

#include <gdal_priv.h>

#include "writer.hpp"

#include "util.hpp"
#include "stats.hpp"

using namespace hlrg::writer;
using namespace geo::util;


GDALWriter::GDALWriter(const std::string& filename, FileType type, int cols, int rows, int bands,
		const std::vector<double>& wavelengths, const std::vector<std::string>& bandNames, char** meta,
		DataType dataType, const std::string& /*interleave*/, const std::string& unit) :
	m_ds(nullptr),
	m_bands(0), m_cols(0), m_rows(0) {

	GDALDataType gtype;
	switch(dataType) {
	case DataType::Byte: gtype = GDT_Byte; break;
	case DataType::Int32: gtype = GDT_Int32; break;
	case DataType::Float32: gtype = GDT_Float32; break;
	default:
		throw std::invalid_argument("Invalid data type.");
	}

	GDALAllRegister();
	CPLSetConfigOption("GDAL_PAM_ENABLED", "NO");

	std::string typeStr = fileTypeAsString(type);
	GDALDriverManager* gm = GetGDALDriverManager();
	GDALDriver* drv = gm->GetDriverByName(typeStr.c_str());
	if(!drv)
		throw std::runtime_error("Driver not found: " + typeStr);

	//CPLStringList options;
	//options.SetNameValue("INTERLEAVE", interleave.c_str());

	m_ds = drv->Create(filename.c_str(), cols, rows, bands, gtype, nullptr);
	if(!m_ds)
		throw std::runtime_error("Failed to open " + filename + " for writing");

	if(meta)
		m_ds->SetMetadata(meta);

	m_bands = m_ds->GetRasterCount();
	m_cols = m_ds->GetRasterXSize();
	m_rows = m_ds->GetRasterYSize();

	if(!bandNames.empty()) {
		for(int i = 1; i <= std::min((int) bandNames.size(), m_bands); ++i) {
			m_ds->GetRasterBand(i)->SetDescription(bandNames[i - 1].c_str());
			m_ds->GetRasterBand(i)->SetMetadataItem("band_name", bandNames[i-1].c_str());
		}
	}
	if(!wavelengths.empty()) {
		m_ds->SetMetadataItem("wavelength_units", unit.c_str());
		for(int i = 1; i <= std::min((int) wavelengths.size(), m_bands); ++i) {
			m_ds->SetMetadataItem("wavelength", std::to_string(wavelengths[i - 1]).c_str());
			m_ds->SetMetadataItem("wavelength_units", unit.c_str());
		}
	}
}

bool GDALWriter::write(const std::vector<double>& buf, int col, int row,
		int cols, int rows, int bufSizeX, int bufSizeY, const std::string& /*id*/) {
	if(col < 0 || col >= m_cols || col + cols > m_cols
			|| row < 0 || row >= m_rows || row + rows > m_rows)
		return false;
	if(bufSizeX <= 0) bufSizeX = cols;
	if(bufSizeY <= 0) bufSizeY = rows;
	const double* data = buf.data();
	for(int i = 1; i <= m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Write, col, row, cols, rows,
				(void*) (data + (i - 1) * bufSizeX * bufSizeY),
				bufSizeX, bufSizeY, GDT_Float64, 0, 0, 0))
			return false;
	}
	return true;
}

void GDALWriter::fill(double v) {
	for(int i = 1; i <= m_bands; ++i)
		m_ds->GetRasterBand(i)->Fill(v);
}

void GDALWriter::fill(int v) {
	fill((double) v);
}

bool GDALWriter::write(const std::vector<int>& buf, int col, int row,
		int cols, int rows, int bufSizeX, int bufSizeY, const std::string& /*id*/) {
	if(col < 0 || col >= m_cols || col + cols > m_cols
			|| row < 0 || row >= m_rows || row + rows > m_rows)
		return false;
	if(bufSizeX <= 0) bufSizeX = cols;
	if(bufSizeY <= 0) bufSizeY = rows;
	const int* data = buf.data();
	for(int i = 1; i <= m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Write, col, row, cols, rows,
				(void*) (data + (i - 1) * bufSizeX * bufSizeY),
				bufSizeX, bufSizeY, GDT_Int32, 0, 0, 0))
			return false;
	}
	return true;
}

void GDALWriter::setProjection(const std::string& projection) {
	m_ds->SetProjection(projection.c_str());
}

void GDALWriter::setTransform(const double* trans) {
	double t[6];
	for(int i = 0; i < 6; ++i)
		t[i] = trans[i];
	m_ds->SetGeoTransform(t);
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
		std::copy_if(buf.begin(), buf.end(), std::back_inserter(values), isnonzero);

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



CSVWriter::CSVWriter(const std::string& filename, const std::vector<double>& wavelengths,
		const std::vector<std::string>& bandNames, const std::string& /*unit*/) {

	m_output.open(filename);
	m_output << "id";

	if(bandNames.empty()) {
		m_output << std::setprecision(3);
		for(double wl : wavelengths)
			m_output << "," << wl;
	} else {
		for(const std::string& bn : bandNames)
			m_output << "," << bn;
	}
	m_output << "\n";
	m_output << std::setprecision(6);
	m_id = 0;
}

std::ofstream& CSVWriter::outstr() {
	return m_output;
}

bool CSVWriter::write(const std::vector<double>& buf, int /*col*/, int /*row*/,
		int /*cols*/, int /*rows*/, int /*bufSizeX*/, int /*bufSizeY*/, const std::string& id) {

	std::string _id(id);
	if(_id.empty())
		_id = std::to_string(++m_id);

	m_output << id;
	for(double v : buf)
		m_output << "," << v;
	m_output << "\n";

	return true;
}

void CSVWriter::fill(double) {

}


void CSVWriter::fill(int) {

}

bool CSVWriter::write(const std::vector<int>& buf, int /*col*/, int /*row*/,
		int /*cols*/, int /*rows*/, int /*bufSizeX*/, int /*bufSizeY*/, const std::string& id) {

	std::string _id(id);
	if(_id.empty())
		_id = std::to_string(++m_id);

	m_output << id;
	for(double v : buf)
		m_output << "," << v;
	m_output << "\n";

	return true;
}


bool CSVWriter::writeStats(const std::string& /*filename*/, const std::vector<std::string>& /*names*/) {
	return true;
}

CSVWriter::~CSVWriter() {
}
