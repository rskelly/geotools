/*
 * reader.hpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#ifndef READER_HPP_
#define READER_HPP_

#include <string>
#include <vector>
#include <map>
#include <unordered_map>

#include <gdal_priv.h>

class Reader {
protected:
	int m_cols;
	int m_rows;
	int m_bands;

	int m_col;
	int m_row;
	int m_bufSize;

	std::map<int, int> m_bandMap;
	int m_minWl; // These are scaled to avoid representation issues.
	int m_maxWl;
	int m_minIdx;
	int m_maxIdx;

public:
	Reader();
	virtual bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows) = 0;
	void setBufSize(int bufSize);
	void setBandMap(const std::map<int, int>& map);
	void setBandRange(double min, double max);
	std::vector<double> getBands() const;
	std::vector<double> getBandRange() const;
	std::vector<int> getIndices() const;
	int bands() const;
	int cols() const;
	int rows() const;
	virtual ~Reader() {}
};

class BandMapReader {
private:
	std::map<int, int> m_bandMap;

public:
	BandMapReader(const std::string& filename, int wlCol, int idxCol, bool hasHeader = true);
	const std::map<int, int>& bandMap() const;
};

class GDALReader : public Reader {
private:
	GDALDataset* m_ds;

public:
	GDALReader(const std::string& filename);
	bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows);
	~GDALReader();
};

class px {
public:
	int c, r;
	std::vector<double> values;
	px(int c, int r) :
		c(c), r(r) {}
	px() : px(0, 0) {}
};

class ROIReader : public Reader {
private:
	std::unordered_map<long, px> m_pixels;

public:
	ROIReader(const std::string& filename);
	bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows);
	~ROIReader();
};

#endif /* READER_HPP_ */
