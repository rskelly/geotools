/*
 * raster.cpp
 *
 *  Created on: Sep 20, 2018
 *      Author: rob
 */

#include "raster.hpp"


using namespace hlrg;

namespace {

	GDALDataType asGDALType(DataType type) {
		switch(type) {
		case UInt16: return GDT_UInt16;
		case Int16: return  GDT_Int16;
		case UInt32: return GDT_UInt32;
		case Int32: return GDT_Int32;
		case Float32: return GDT_Float32;
		case Float64: return GDT_Float64;
		case Byte: return GDT_Byte;
		}
		throw std::runtime_error("Unknown type.");
	}

	template <class T>
	GDALDataType asGDALType(T v) {
		size_t t = typeid(T).hash_code();
		if(t == typeid(char).hash_code()) {
			return GDT_Byte;
		} else if(t == typeid(int).hash_code()) {
			return GDT_Int32;
		} else if(t == typeid(unsigned int).hash_code()) {
			return GDT_UInt32;
		} else if(t == typeid(float).hash_code()) {
			return GDT_Float32;
		} else if(t == typeid(double).hash_code()) {
			return GDT_Float64;
		} else if(t == typeid(short).hash_code()) {
			return GDT_Int16;
		} else if(t == typeid(unsigned short).hash_code()) {
			return GDT_UInt16;
		}
		throw std::runtime_error("Unknown type.");
	}

}

Raster::Raster(const std::string& filename) :
	m_ds(nullptr),
	m_bands(0),
	m_cols(0),
	m_rows(0),
	m_row(-1) {

	GDALAllRegister();

	if(!(m_ds = (GDALDataset*) GDALOpen(filename.c_str(), GA_ReadOnly)))
		throw std::invalid_argument("Failed to open dataset.");

	m_bands = m_ds->GetRasterCount();
	m_cols = m_ds->GetRasterXSize();
	m_rows = m_ds->GetRasterYSize();
}

Raster::Raster(const std::string& filename, int cols, int rows, int bands, int srid, DataType type, Raster* parent) :
	m_ds(nullptr),
	m_bands(0),
	m_cols(0),
	m_rows(0),
	m_row(-1) {

	GDALAllRegister();

	static char op[] {"INTERLEAVE=BSQ"};
	static char* op0 = &op[0];

	GDALDriverManager* dm = GetGDALDriverManager();
	GDALDriver* drv = dm->GetDriverByName("ENVI");
	GDALDataType gtype = asGDALType(type);
	if(!(m_ds = (GDALDataset*) drv->Create(filename.c_str(), cols, rows, bands, gtype, &op0)))
		throw std::invalid_argument("Failed to create dataset.");

	m_bands = m_ds->GetRasterCount();
	m_cols = m_ds->GetRasterXSize();
	m_rows = m_ds->GetRasterYSize();

	if(parent) {
		char** meta = parent->metadata();
		m_ds->SetMetadata(meta);
		for(int b = 1; b <= m_bands; ++b) {
			char** bmeta = parent->band(b)->GetMetadata();
			m_ds->GetRasterBand(b)->SetMetadata(bmeta);
		}
	}
}

char** Raster::metadata() const {
	return m_ds->GetMetadata();
}

GDALDataset* Raster::ds() const {
	return m_ds;
}

GDALRasterBand* Raster::band(int b) const {
	return m_ds->GetRasterBand(b);
}

template <class T>
bool Raster::next(std::vector<T>& buf) {
	T tmp = 0;
	GDALDataType gtype = asGDALType(tmp);
	buf.resize(m_bands * m_cols);
	std::fill(buf.begin(), buf.end(), 0);
	if(++m_row < m_rows) {
		for(int i = 1; i <= m_bands; ++i) {
			GDALRasterBand* band = m_ds->GetRasterBand(i);
			if(CE_None != band->RasterIO(GF_Read, 0, m_row, m_cols, 1, (T*) (buf.data() + (i - 1) * m_cols), m_cols, 1, gtype, 0, 0, nullptr))
				throw std::runtime_error("Failed to read from raster.");
		}
		return true;
	}
	return false;
}

template <class T>
bool Raster::write(std::vector<T>& buf, int row) {
	T tmp = 0;
	GDALDataType gtype = asGDALType(tmp);
	if(row < 0 || row >= m_rows)
		return false;
	for(int i = 1; i <= m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(CE_None != band->RasterIO(GF_Write, 0, row, m_cols, 1, (T*) (buf.data() + (i - 1) * m_cols), m_cols, 1, gtype, 0, 0, nullptr))
			throw std::runtime_error("Failed to write to raster.");
	}
	return true;
}

template <class T>
bool Raster::get(std::vector<T>& buf, int row) {
	T tmp = 0;
	GDALDataType gtype = asGDALType(tmp);
	buf.resize(m_bands * m_cols);
	std::fill(buf.begin(), buf.end(), 0);
	if(row < 0 || row >= m_rows)
		return false;
	for(int i = 1; i <= m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(CE_None != band->RasterIO(GF_Read, 0, row, m_cols, 1, (T*) (buf.data() + (i - 1) * m_cols), m_cols, 1, gtype, 0, 0, nullptr))
			return false;
	}
	return true;
}
int Raster::bands() const {
	return m_bands;
}

int Raster::cols() const {
	return m_cols;
}

int Raster::rows() const {
	return m_rows;
}

void Raster::reset() {
	m_row = -1;
}

Raster::~Raster() {
	GDALClose(m_ds);
}


template bool Raster::get<int>(std::vector<int>& buf, int row);
template bool Raster::get<unsigned int>(std::vector<unsigned int>& buf, int row);
template bool Raster::get<short>(std::vector<short>& buf, int row);
template bool Raster::get<unsigned short>(std::vector<unsigned short>& buf, int row);
template bool Raster::get<char>(std::vector<char>& buf, int row);
template bool Raster::get<float>(std::vector<float>& buf, int row);
template bool Raster::get<double>(std::vector<double>& buf, int row);

template bool Raster::write<int>(std::vector<int>& buf, int row);
template bool Raster::write<unsigned int>(std::vector<unsigned int>& buf, int row);
template bool Raster::write<short>(std::vector<short>& buf, int row);
template bool Raster::write<unsigned short>(std::vector<unsigned short>& buf, int row);
template bool Raster::write<char>(std::vector<char>& buf, int row);
template bool Raster::write<float>(std::vector<float>& buf, int row);
template bool Raster::write<double>(std::vector<double>& buf, int row);

template bool Raster::next<int>(std::vector<int>& buf);
template bool Raster::next<unsigned int>(std::vector<unsigned int>& buf);
template bool Raster::next<short>(std::vector<short>& buf);
template bool Raster::next<unsigned short>(std::vector<unsigned short>& buf);
template bool Raster::next<char>(std::vector<char>& buf);
template bool Raster::next<float>(std::vector<float>& buf);
template bool Raster::next<double>(std::vector<double>& buf);


