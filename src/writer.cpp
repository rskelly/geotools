/*
 * writer.cpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#include <gdal_priv.h>

#include "writer.hpp"

GDALWriter::GDALWriter(const std::string& filename, int cols, int rows, int bands) :
	m_ds(nullptr),
	m_bands(0), m_cols(0), m_rows(0) {

	GDALAllRegister();
	GDALDriverManager gm;
	GDALDriver* drv = gm.GetDriverByName(filename.c_str());
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
	for(int i = 1; i < m_bands; ++i) {
		GDALRasterBand* band = m_ds->GetRasterBand(i);
		if(band->RasterIO(GF_Write, col, row, cols, rows,
				(void*) (data + (i - 1) * bufSize * bufSize),
				bufSize, bufSize, GDT_Float64, 0, 0, 0))
			return false;
	}
	return true;
}

GDALWriter::~GDALWriter() {
	GDALClose(m_ds);
}

