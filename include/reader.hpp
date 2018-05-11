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
public:
	virtual bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows) = 0;
	virtual void setBandMap(std::map<int, int>& map) = 0;
	virtual void setBandMap(std::string& bandfile) = 0;
	virtual void setBandRange(double min, double max) = 0;
	virtual std::vector<double> getBands() const = 0;
	virtual std::vector<double> getBandRange() const = 0;
	virtual std::vector<int> getIndices() const = 0;
	virtual int bands() const = 0;
	virtual int cols() const = 0;
	virtual int rows() const = 0;
	virtual ~Reader() {}
};

class GDALReader : public Reader {
private:
	GDALDataset* m_ds;
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
	GDALReader(const std::string& filename);
	void setBufSize(int bufSize);
	bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows);
	void setBandMap(std::map<int, int>& map);
	void setBandMap(std::string& bandfile);
	void setBandRange(double min, double max);
	std::vector<double> getBands() const;
	std::vector<double> getBandRange() const;
	std::vector<int> getIndices() const;
	int bands() const;
	int cols() const;
	int rows() const;
	~GDALReader();
};

class px {
	int c, r;
	std::vector<double> values;
	px(int c, int r) :
		c(c), r(r) {}
	px() : px(0, 0) {}
};

class TextReader : public Reader {
private:
	int m_cols;
	int m_rows;
	int m_bands;
	std::unordered_map<long, px> m_pixels;
	int m_col;
	int m_row;
	int m_bufSize;

	std::map<int, int> m_bandMap;
	int m_minWl; // These are scaled to avoid representation issues.
	int m_maxWl;
	int m_minIdx;
	int m_maxIdx;

public:
	TextReader(const std::string& filename);
	void setBufSize(int bufSize);
	bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows);
	void setBandMap(std::map<int, int>& map);
	void setBandMap(std::string& bandfile);
	void setBandRange(double min, double max);
	std::vector<double> getBands() const;
	std::vector<double> getBandRange() const;
	std::vector<int> getIndices() const;
	int bands() const;
	int cols() const;
	int rows() const;
	~TextReader();
};

#endif /* READER_HPP_ */
