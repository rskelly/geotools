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

#include <gdal_priv.h>

class Reader {
public:
	virtual bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows) = 0;
	virtual void setBandMap(std::map<double, int>& map) = 0;
	virtual void setBandRange(double min, double max) = 0;
	virtual std::vector<double> getBands() const = 0;
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

	std::map<double, int> m_bandMap;
	double m_minWl;
	double m_maxWl;
	int m_minIdx;
	int m_maxIdx;

public:
	GDALReader(const std::string& filename);
	void setBufSize(int bufSize);
	bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows);
	void setBandMap(std::map<double, int>& map);
	void setBandRange(double min, double max);
	std::vector<double> getBands() const;
	std::vector<int> getIndices() const;
	int bands() const;
	int cols() const;
	int rows() const;
	~GDALReader();
};
#endif /* READER_HPP_ */
