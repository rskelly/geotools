/*
 * writer.hpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#ifndef WRITER_HPP_
#define WRITER_HPP_

#include <vector>
#include <string>

#include <gdal_priv.h>

class Writer {
public:
	virtual bool write(std::vector<double>& buf, int col, int row, int cols, int rows, int bufSize) = 0;
	virtual bool writeStats(const std::string& filename) = 0;
	virtual ~Writer() {}
};

class GDALWriter {
private:
	GDALDataset* m_ds;
	int m_bands;
	int m_cols;
	int m_rows;

public:
	GDALWriter(const std::string& filename, int cols, int rows, int bands, const std::string& fieldName = "", const std::vector<std::string>& bandNames = {});
	bool write(std::vector<double>& buf, int col, int row, int cols, int rows, int bufSize);
	bool writeStats(const std::string& filename);
	~GDALWriter();
};


#endif /* WRITER_HPP_ */
