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

namespace hlrg {

/**
 * Abstract class for object that read spectral datasets.
 */
class Reader {
protected:
	int m_cols;		///<! Number of columns.
	int m_rows; 	///<! Number of rows.
	int m_bands;	///<! Number of bands.

	int m_col;		///<! Current column for iteration.
	int m_row;		///<! Current row for iteration.
	int m_bufSize;	///<! The buffer size.

	std::map<int, int> m_bandMap;			///<! A index-to-band mapping.
	std::vector<std::string> m_bandNames;	///<! A list of band names; possibly for storing a representation of the band wavelength, etc.
	int m_minWl; 							///<! Minimum wavelength; these are scaled to avoid representation issues.
	int m_maxWl;							///<! Maximum wavelength.
	int m_minIdx;							///<! Minimum band index.
	int m_maxIdx;							///<! Maximum band index.

public:
	/**
	 * Default constructor.
	 */
	Reader();

	/**
	 * Fill the buffer with the next available row of data. Data will be stored sequentially by band,
	 * so band 1 will be stored as pixel 0-n, and then band 2, pixel 0-n, etc. The column
	 * and row references are updated with the current values representing the region just read.
	 *
	 * @param buf A buffer large enough to store the pixels and bands of a row, in band-sequential order.
	 * @param col A reference to a value that will be updated with the current column.
	 * @param col A reference to a value that will be updated with the current row.
	 * @param col A reference to a value that will be updated with the number of columns in the buffer.
	 * @param col A reference to a value that will be updated with the numbe of rows in the buffer.
	 * @return True if there is another row to be read after this one.
	 */
	virtual bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows) = 0;

	/**
	 * Set the size of the buffer for reading.
	 *
	 * @param bufSize The size of the buffer.
	 */
	void setBufSize(int bufSize);

	/**
	 * Set the band map; a mapping of band index to band number.
	 *
	 * @param map A map containing the mapping from band index to band number.
	 */
	void setBandMap(const std::map<int, int>& map);

	/**
	 * Sets the range of wavelengths.
	 *
	 * @param min The minimum wavelength.
	 * @param max The maximum wavelength.
	 */
	void setBandRange(double min, double max);

	/**
	 * Return a vector containing the wavelengths.
	 *
	 * @return A vector containing the wavelengths.
	 */
	std::vector<double> getWavelengths() const;

	/**
	 * Returns a vector containing the band names.
	 *
	 * @return A vector containing the band names.
	 */
	std::vector<std::string> getBandNames() const;

	/**
	 * Returns a two-element vector containing the min and max wavelengths.
	 *
	 * @return A two-element vector containing the min and max wavelengths.
	 */
	std::vector<double> getBandRange() const;

	/**
	 * Returns a two-element vector containing the min and max indices.
	 *
	 * @return A two-element vector containing the min and max indices.
	 */
	std::vector<int> getIndices() const;

	/**
	 * Returns the number of bands.
	 *
	 * @return The number of bands.
	 */
	int bands() const;

	/**
	 * Returns the number of columns.
	 *
	 * @return The number of columns.
	 */
	int cols() const;

	/**
	 * Returns the number of rows.
	 *
	 * @return The number of rows.
	 */
	int rows() const;

	virtual ~Reader() {}
};


/**
 * Reads a file containing a band map and provides the result as a map.
 */
class BandMapReader {
private:
	std::map<int, int> m_bandMap;

public:
	/**
	 * Construct the band map.
	 *
	 * @param filename The source of the band map data.
	 * @param wlCol The column containing wavelengths.
	 * @param idxCol The column containing the band indices.
	 * @param hasHeader True if the source has a header row which should be ignored.
	 */
	BandMapReader(const std::string& filename, int wlCol, int idxCol, bool hasHeader = true);

	/**
	 * Return a reference to the internal band map.
	 *
	 * @return The band map.
	 */
	const std::map<int, int>& bandMap() const;
};

/**
 * An implementation of Reader that can read from GDAL data sources.
 */
class GDALReader : public Reader {
private:
	GDALDataset* m_ds;

public:
	/**
	 * Construct the reader around the given filename.
	 *
	 * @param filename The filename of a GDAL data source.
	 */
	GDALReader(const std::string& filename);

	/**
	 * Fill the buffer with the next available row of data. Data will be stored sequentially by band,
	 * so band 1 will be stored as pixel 0-n, and then band 2, pixel 0-n, etc. The column
	 * and row references are updated with the current values representing the region just read.
	 *
	 * @param buf A buffer large enough to store the pixels and bands of a row, in band-sequential order.
	 * @param col A reference to a value that will be updated with the current column.
	 * @param col A reference to a value that will be updated with the current row.
	 * @param col A reference to a value that will be updated with the number of columns in the buffer.
	 * @param col A reference to a value that will be updated with the numbe of rows in the buffer.
	 * @return True if there is another row to be read after this one.
	 */
	bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows);

	 ~GDALReader();
};

/**
 * A utility class representing a pixel.
 */
class px {
public:
	int c, r;
	std::vector<double> values;
	px(int c, int r) :
		c(c), r(r) {}
	px() : px(0, 0) {}
};


/**
 * An implementation of Reader that reads ENVI ROI files.
 */
class ROIReader : public Reader {
private:
	std::unordered_map<long, px> m_pixels;

public:
	/**
	 * Construct the reader around the given ROI file.
	 *
	 * @param filename An ROI file.
	 */
	ROIReader(const std::string& filename);

	/**
	 * Fill the buffer with the next available row of data. Data will be stored sequentially by band,
	 * so band 1 will be stored as pixel 0-n, and then band 2, pixel 0-n, etc. The column
	 * and row references are updated with the current values representing the region just read.
	 *
	 * @param buf A buffer large enough to store the pixels and bands of a row, in band-sequential order.
	 * @param col A reference to a value that will be updated with the current column.
	 * @param col A reference to a value that will be updated with the current row.
	 * @param col A reference to a value that will be updated with the number of columns in the buffer.
	 * @param col A reference to a value that will be updated with the numbe of rows in the buffer.
	 * @return True if there is another row to be read after this one.
	 */
	bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows);

	~ROIReader();
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

	/**
	 * Get the frame corresponding to the given timestamp.
	 * Returns true if the frame was found.
	 *
	 * @param utcTime The timestamp to search for.
	 * @param frame A value that will be updated with the index of the frame.
	 * @return True if the frame is found, false otherwise.
	 */
	bool getFrame(long utcTime, int& frame) const;

	/**
	 * Get the frame nearest to the given timestamp.
	 * Returns true if a frame was found.
	 *
	 * @param utcTime The timestamp to search for.
	 * @param actualUtcTime A value that will be updated with the nearest timestamp to the one given.
	 * @param frame A value that will be updated with the index of the frame.
	 * @return True if a frame is found, false otherwise.
	 */
	bool getNearestFrame(long utcTime, long& actualUtcTime, int& frame) const;

	/**
	 * Get the timestamp corresponding to the given frame.
	 * Returns true if the timestamp was found.
	 *
	 * @param frame The frame to search for.
	 * @param utcTime A value that will be updated with the timestamp.
	 * @return True if the timestamp is found, false otherwise.
	 */
	bool getTime(int frame, long& utcTime) const;

	/**
	 * Get the timestamp nearest to the given frame.
	 * Returns true if a timestamp was found.
	 *
	 * @param frame The frame to search for.
	 * @param actualFrame A value that will be updated with the nearest frame to the one given.
	 * @param utcTime A value that will be updated with the index of the timestamp.
	 * @return True if a timestmap is found, false otherwise.
	 */
	bool getNearestTime(int frame, int& actualFrame, long& utcTime) const;

	~FrameIndexReader();
};

/**
 * Represents a row in the imu_gps file from the APX-15 IMU.
 */
class IMUGPSRow {
public:
	double roll;	///<! The platform roll.
	double pitch;	///<! The platform pitch.
	double yaw;		///<! The platform yaw.
	double lat;		///<! The platform latitude.
	double lon;		///<! The platform longitude.
	double alt;		///<! The platform altitude.
	long gpsTime;	///<! The GPS timestamp.
	long utcTime;	///<! The UTC timestamp (since epoch, ms).
	int status;		///<! The status.
	double heading;	///<! The platform heading.

	size_t index;	///<! An index for this row in relation to other rows.

	/***
	 * Construct an IMUGPSRow by passing an input stream from which it can read information.
	 * The offset argument is for adjusting the timestamps according to a known offset, in
	 * milliseconds.
	 *
	 * @param in An input stream.
	 * @param msOffset A time offset to apply to the times in the rows, in milliseconds.
	 */
	IMUGPSRow(std::istream& in, double msOffset);
};

/**
 * Loads the APX-15 IMU/GPS table from a text file.
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
	 * @param gpsTime The timestamp as emitted by the APX.
	 * @param utcTime A value that will be update with the UTC timestamp.
	 * @return True if the utcTime value has been set successfully.
	 */
	bool getUTCTime(long gpsTime, long& utcTime);

	/**
	 * Compute and return the interpolated GPS timestamp corresponding to the given
	 * UTC timestamp.
	 *
	 * @param utcTime The UTC timestamp.
	 * @param gpsTime A value that will be updated with the GPS timestamp.
	 * @return True if the gpsTime is updated successfully.
	 */
	bool getGPSTime(long utcTime, long& gpsTime);

	~IMUGPSReader();

};


/**
 * Represents a row from the OceanOptics Flame data output.
 */
class FlameRow {
public:
	long dateTime;						///<! A timestamp (ms) derived from the date field.
	long utcTime;						///<! The UTC timestamp (ms) as output.
	std::vector<double> bands;			///<! A vector containing the band values.
	std::vector<double> wavelengths;	///<! A vector containing the wavelengths. TODO: This should be stored as ints.

	/**
	 * Read the row data from the given input stream.
	 * The offset argument applies a known offst (in ms) to the times in the row.
	 *
	 * @param in An input stream.
	 * @param msOffset A time offset in milliseconds.
	 * @return True if the row was read successfully.
	 */
	bool read(std::istream& in, double msOffset);
};

/**
 * Reads the CSV for a convolved Flame dataset, convolved using the convolve program
 * in this toolset. The row format is:
 *
 * 	date,timestamp,[band wl1],[band wl2],[etc...]
 */
class FlameReader {
private:
	std::ifstream m_in;		///<! The file reader.
	double m_msOffset;		///<! The time offset in milliseconds.
	std::string m_filename;	///<! The data file name.

public:
	std::vector<double> wavelengths;	///<! The list of wavelengths. TODO: This should be stored as ints.

	/**
	 * Construct a FlameReader using the given filename and time offset.
	 *
	 * @param filename The filename of the Flame output dataset.
	 * @param msOffset A time offset in milliseconds to apply to the times stored in each row.
	 */
	FlameReader(const std::string& filename, double msOffset);

	/**
	 * Return the number of rows in the file. This requires reading through the whole file.
	 *
	 * @return The number of rows in the file.
	 */
	int rows();

	/**
	 * Read the next row of data into the given row object.
	 *
	 * @param row A FlameRow instance to populate with values from the next row.
	 * @return True if there is another row to read after the present row.
	 */
	bool next(FlameRow& row);
};

} // hlrg

#endif /* READER_HPP_ */
