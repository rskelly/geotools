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
#include <fstream>

#include <gdal_priv.h>

#include "bintree.hpp"

#define MIN_VALUE 0.000001
#define WL_SCALE 100000

class Reader {
protected:
	int m_cols;
	int m_rows;
	int m_bands;

	int m_col;
	int m_row;
	int m_bufSize;

	std::map<int, int> m_bandMap;
	std::vector<std::string> m_bandNames;
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
	std::vector<double> getWavelengths() const;
	std::vector<std::string> getBandNames() const;
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


/**
 * Provides methods for reading rows out of a raster.
 */
class Raster {
private:
	GDALDataset* m_ds;
	int m_bands;
	int m_cols;
	int m_rows;
	int m_row;

public:
	/**
	 * Construct a raster object using the given filename.
	 */
	Raster(const std::string& filename);

	int bands() const;

	int cols() const;

	int rows() const;

	/**
	 * Read all bands of the next row. This will write the pixel rows for each band
	 * end-to-end into the given buffer.
	 *
	 * @param buf A vector which contains the pixel row for each band, end to end.
	 * @return True if there is another row to be read. False otherwise.
	 */
	bool next(std::vector<double>& buf);

	/**
	 * Get the row by zero-based index.
	 *
	 * @param buf A vector which contains the pixel row for each band, end to end.
	 * @param row The zero-based row index.
	 * @return True if the read was successful.
	 */
	bool get(std::vector<uint16_t>& buf, int row);

	/**
	 * Reset the file pointer to the beginning.
	 */
	void reset();

	~Raster();
};


/**
 * Reads the frame index time files from the Hyperspect Nano instrument.
 */
class FrameIndexReader {
private:
	BinTree<long, int>* m_frames;
	BinTree<int, long>* m_times;
public:

	/**
	 * Loads the frame index into a map. Indexed by frame index (0-based, value is the GPS timestamp in us).
	 *
	 * @param filename The filename of the frame index file.
	 */
	FrameIndexReader(const std::string& filename);

	bool getFrame(long time, int& frame) const;

	bool getNearestFrame(long time, long& actualTime, int& frame) const;

	bool getTime(int frame, long& time) const;

	bool getNearestTime(int frame, int& actualFrame, long& time) const;

	~FrameIndexReader();
};

class IMUGPSRow {
public:
	double roll;
	double pitch;
	double yaw;
	double lat;
	double lon;
	double alt;
	long timestamp;
	long date; // As a utc timestamp
	int status;
	double heading;

	size_t index;

	IMUGPSRow(std::ifstream& in, double msOffset);
};

/**
 * Loads the IMU/GPS table from a text file.
 */
class IMUGPSReader {
private:
	std::ifstream m_in;
	BinTree<long, IMUGPSRow*>* m_gpsTimesT;
	BinTree<long, IMUGPSRow*>* m_utcTimesT;
	std::vector<IMUGPSRow*> m_rows;
	size_t m_lastIndex;

public:

	/**
	 * Load the file.
	 *
	 * @param filename The filename.
	 */
	IMUGPSReader(const std::string& filename, double msOffset);

	/**
	 * Compute and return the interpolated UTC timestamp since the epoch (Jan 1, 1970) in miliseconds.
	 *
	 * @param timestamp The timestamp as emitted by the APX.
	 * @return The UTC timestamp in miliseconds.
	 */
	bool getUTCTime(long timestamp, long& utcTime);

	bool getGPSTime(long date, long& gpsTime);

	~IMUGPSReader();

};

class FlameRow {
public:
	long date;
	long timestamp;
	std::vector<double> bands;
	std::vector<double> wavelengths;

	bool read(std::ifstream& in, double msOffset);
};

/**
 * Reads the CSV for a convolved Flame dataset. The row format is,
 * date,timestamp,[band wl1],[band wl2],[etc...]
 */
class FlameReader {
private:
	std::ifstream m_in;
	double m_msOffset;
public:
	std::vector<double> wavelengths;

	FlameReader(const std::string& filename, double msOffset);

	bool next(FlameRow& row);
};

#endif /* READER_HPP_ */
