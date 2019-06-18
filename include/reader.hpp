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
#include "util.hpp"
#include "ds/kdtree.hpp"

namespace hlrg {
namespace reader {

constexpr double MIN_VALUE = 0.000001; // Note: Can't use this; screws up hull std::numeric_limits<double>::min();
constexpr double WL_SCALE = 100000;

/**
 * Abstract class for object that read spectral datasets.
 */
class Reader {
protected:
	std::string m_filename;
	int m_cols;								///<! Number of columns.
	int m_rows; 							///<! Number of rows.
	int m_bands;							///<! Number of bands.

	int m_col;								///<! Current column for iteration.
	int m_row;								///<! Current row for iteration.
	int m_bufSize;							///<! The buffer size.

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
	 * Return the FileType represented by this reader.
	 *
	 * \return The FileType represented by this reader.
	 */
	hlrg::util::FileType fileType() const;

	/**
	 * Fill the buffer with the next available row of data. Data will be stored sequentially by band,
	 * so band 1 will be stored as pixel 0-n, and then band 2, pixel 0-n, etc. The column
	 * and row references are updated with the current values representing the region just read.
	 *
	 * \param[out] id A string to hold the identifier for the item, if there is one.
	 * \param[out] buf A buffer to store the pixels and bands of a row, in band-sequential order.
	 * \param[out] cols A reference to a value that will be updated with the number of columns in the buffer.
	 * \param[out] col A reference to the column index that was read.
	 * \param[out] row A reference to the row index that was read.
	 * \return True if the row was read successfully.
	 */
	virtual bool next(std::string& id, std::vector<double>& buf, int& cols, int& col, int& row) = 0;

	/**
	 * Set the size of the buffer for reading.
	 *
	 * \param bufSize The size of the buffer.
	 */
	void setBufSize(int bufSize);

	/**
	 * Set the band map; a mapping of band index to band number.
	 *
	 * \param map A map containing the mapping from band index to band number.
	 */
	void setBandMap(const std::map<int, int>& map);

	/**
	 * Sets the range of wavelengths.
	 *
	 * \param min The minimum wavelength.
	 * \param max The maximum wavelength.
	 */
	void setBandRange(double min, double max);

	/**
	 * Return a vector containing the wavelengths.
	 *
	 * \return A vector containing the wavelengths.
	 */
	std::vector<double> getWavelengths() const;

	/**
	 * Returns a vector containing the band names.
	 *
	 * \return A vector containing the band names.
	 */
	std::vector<std::string> getBandNames() const;

	/**
	 * Returns a two-element vector containing the min and max wavelengths.
	 *
	 * \return A two-element vector containing the min and max wavelengths.
	 */
	std::vector<double> getBandRange() const;

	/**
	 * Returns a two-element vector containing the min and max indices.
	 *
	 * \return A two-element vector containing the min and max indices.
	 */
	std::vector<int> getIndices() const;

	/**
	 * Returns the mapping of band number (1-based) to wavelength. Wavelength
	 * is a floating point number scaled by WL_SCALE to be an integer.
	 *
	 * \return A mapping of band number to wavelength, scaled to an integer.
	 */
	std::map<int, int> getBandMap();

	/**
	 * Returns the number of bands.
	 *
	 * \return The number of bands.
	 */
	int bands() const;

	/**
	 * Returns the number of columns.
	 *
	 * \return The number of columns.
	 */
	int cols() const;

	/**
	 * Returns the number of rows.
	 *
	 * \return The number of rows.
	 */
	int rows() const;

	virtual ~Reader() {}
};

class Point {
private:
	double m_x, m_y, m_z;
	std::string m_id;

public:
	Point(double x = 0, double y = 0, double z = 0) :
		m_x(x), m_y(y), m_z(z) {}

	double operator[](int idx) const {
		switch(idx % 3) {
		case 0: return m_x;
		case 1: return m_y;
		default: return m_z;
		}
	}

	double& operator[](int idx) {
		switch(idx % 3) {
		case 0: return m_x;
		case 1: return m_y;
		default: return m_z;
		}
	}

	double x() const { return m_x; }
	double y() const { return m_y; }
	double z() const { return m_z; }

};


/**
 * Reads and stores a set of points from a spatial database and stores them
 * in a tree where they can be searched by radius.
 */
class PointSetReader {
private:
	geo::ds::KDTree<hlrg::reader::Point>* m_tree;

public:
	PointSetReader(const std::string& filename, const std::string& layer);
	static std::vector<std::string> getLayerNames(const std::string& filename);
	int search(double x, double y, double radius, std::vector<hlrg::reader::Point*>& pts);
	~PointSetReader();
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
	 * \param filename The source of the band map data.
	 * \param wlCol The column containing wavelengths.
	 * \param idxCol The column containing the band indices.
	 * \param hasHeader True if the source has a header row which should be ignored.
	 */
	BandMapReader(const std::string& filename, int wlCol, int idxCol, bool hasHeader = true);

	/**
	 * Return a reference to the internal band map.
	 *
	 * \return The band map.
	 */
	const std::map<int, int>& bandMap() const;
};

/**
 * An implementation of Reader that can read from GDAL data sources.
 */
class GDALReader : public Reader {
private:
	GDALDataset* m_ds;

	size_t m_mappedSize;		///<! The size of mapped memory for remapping an interleaved raster to a list of spectra.
	size_t m_mappedMinBand;		///<! The first mapped band (1-based).
	double* m_mapped;			///<! The pointer to mapped memory for remapping an interleaved raster to a list of spectra.
	int m_mappedBands;			///<! The number of bands mapped into memory.
	std::unique_ptr<hlrg::util::TmpFile> m_mappedFile;
	double m_trans[6];

	void loadBandMap();

public:
	/**
	 * Construct the reader around the given filename.
	 *
	 * \param filename The filename of a GDAL data source.
	 */
	GDALReader(const std::string& filename);

	double toX(int col) const;
	double toY(int row) const;
	double toCol(double x) const;
	double toRow(double y) const;
	double resX() const;
	double resY() const;

	/**
	 * Remap the raster into a memory-mapped list of spectra, organized by
	 * row/col/band. This is freed when the reader is destroyed.
	 *
	 * \param minWl the minimum wavelength of the mapped region.
	 * \param maxWl the maximum wavelength of the mapped region.
	 */
	void remap(double minWl, double maxWl);

	/**
	 * Remap the raster into a memory-mapped list of spectra, organized by
	 * row/col/band. This is freed when the reader is destroyed.
	 *
	 * \param minBand the first band of the mapped region.
	 * \param maxBand the last band of the mapped region.
	 */
	void remap(int minBand, int maxBand);

	/**
	 * Remap the raster into a memory-mapped list of spectra, organized by
	 * row/col/band. This is freed when the reader is destroyed.
	 */
	void remap();

	/**
	 * Get the value of the mapped spectrum at the column and row for the given wavelength.
	 *
	 * \param col The column.
	 * \param row The row.
	 * \param wl The wavelength.
	 * \return The intensity value. NaN if out of range.
	 */
	double mapped(int col, int row, double wl);

	/**
	 * Get the value of the mapped spectrum at the column and row for the given wavelength.
	 *
	 * \param col The column.
	 * \param row The row.
	 * \param The band (1-based).
	 * \return The intensity value. NaN if out of range.
	 */
	double mapped(int col, int row, int band);

	/**
	 * Get the entire mapped spectrum at the given column and row.
	 *
	 * \param col The column.
	 * \param row The row.
	 * \param values A vector to contain the spectrum.
	 * \return True if successful.
	 */
	bool mapped(int col, int row, std::vector<double>& values);

	int toCol(double x);

	int toRow(double y);

	int getInt(double x, double y);

	int getInt(int col, int row);

	float getFloat(double x, double y);

	float getFloat(int col, int row);

	bool next(std::string& id, std::vector<double>& buf, int& cols, int& col, int& row);

	bool next(std::vector<double>& buf, int band, int& cols, int& col, int& row);

	 ~GDALReader();
};

namespace {
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

}

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
	 * \param filename An ROI file.
	 */
	ROIReader(const std::string& filename);

	/**
	 * Fill the buffer with the next available row of data. Data will be stored sequentially by band,
	 * so band 1 will be stored as pixel 0-n, and then band 2, pixel 0-n, etc. The column
	 * and row references are updated with the current values representing the region just read.
	 *
	 * \param buf A buffer large enough to store the pixels and bands of a row, in band-sequential order.
	 * \param col A reference to a value that will be updated with the current column.
	 * \param col A reference to a value that will be updated with the current row.
	 * \param col A reference to a value that will be updated with the number of columns in the buffer.
	 * \param col A reference to a value that will be updated with the numbe of rows in the buffer.
	 * \return True if the row was read successfully.
	 */
	bool next(std::vector<double>& buf, int& col, int& row, int& cols, int& rows);

	~ROIReader();
};


/**
 * Reads the frame index time files from the Hyperspect Nano instrument.
 */
class FrameIndexReader {
private:
	hlrg::ds::BinTree<long, int>* m_frames;
	hlrg::ds::BinTree<int, long>* m_times;
public:

	/**
	 * Loads the frame index into a map. Indexed by frame index (0-based, value is the GPS timestamp in us).
	 *
	 * \param filename The filename of the frame index file.
	 */
	FrameIndexReader(const std::string& filename);

	/**
	 * Get the frame corresponding to the given timestamp.
	 * Returns true if the frame was found.
	 *
	 * \param utcTime The timestamp to search for.
	 * \param frame A value that will be updated with the index of the frame.
	 * \return True if the frame is found, false otherwise.
	 */
	bool getFrame(long utcTime, int& frame) const;

	/**
	 * Get the frame nearest to the given timestamp.
	 * Returns true if a frame was found.
	 *
	 * \param utcTime The timestamp to search for.
	 * \param actualUtcTime A value that will be updated with the nearest timestamp to the one given.
	 * \param frame A value that will be updated with the index of the frame.
	 * \return True if a frame is found, false otherwise.
	 */
	bool getNearestFrame(long utcTime, long& actualUtcTime, int& frame) const;

	/**
	 * Get the timestamp corresponding to the given frame.
	 * Returns true if the timestamp was found.
	 *
	 * \param frame The frame to search for.
	 * \param utcTime A value that will be updated with the timestamp.
	 * \return True if the timestamp is found, false otherwise.
	 */
	bool getTime(int frame, long& utcTime) const;

	/**
	 * Get the timestamp nearest to the given frame.
	 * Returns true if a timestamp was found.
	 *
	 * \param frame The frame to search for.
	 * \param actualFrame A value that will be updated with the nearest frame to the one given.
	 * \param utcTime A value that will be updated with the index of the timestamp.
	 * \return True if a timestmap is found, false otherwise.
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
	 * \param in An input stream.
	 * \param msOffset A time offset to apply to the times in the rows, in milliseconds.
	 */
	IMUGPSRow(std::istream& in, double msOffset);
};

/**
 * Loads the APX-15 IMU/GPS table from a text file.
 */
class IMUGPSReader {
private:
	std::ifstream m_in;
	hlrg::ds::BinTree<long, IMUGPSRow*>* m_gpsTimesT;
	hlrg::ds::BinTree<long, IMUGPSRow*>* m_utcTimesT;
	std::vector<IMUGPSRow*> m_rows;
	size_t m_lastIndex;

public:

	/**
	 * Load the file.
	 *
	 * \param filename The filename.
	 */
	IMUGPSReader(const std::string& filename, double msOffset);

	/**
	 * Compute and return the interpolated UTC timestamp since the epoch (Jan 1, 1970) in miliseconds.
	 *
	 * \param gpsTime The timestamp as emitted by the APX.
	 * \param utcTime A value that will be update with the UTC timestamp.
	 * \return True if the utcTime value has been set successfully.
	 */
	bool getUTCTime(long gpsTime, long& utcTime);

	/**
	 * Compute and return the interpolated GPS timestamp corresponding to the given
	 * UTC timestamp.
	 *
	 * \param utcTime The UTC timestamp.
	 * \param gpsTime A value that will be updated with the GPS timestamp.
	 * \return True if the gpsTime is updated successfully.
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
	 * \param in An input stream.
	 * \param msOffset A time offset in milliseconds.
	 * \return True if the row was read successfully.
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
	 * \param filename The filename of the Flame output dataset.
	 * \param msOffset A time offset in milliseconds to apply to the times stored in each row.
	 */
	FlameReader(const std::string& filename, double msOffset);

	/**
	 * Return the number of rows in the file. This requires reading through the whole file.
	 *
	 * \return The number of rows in the file.
	 */
	int rows();

	/**
	 * Read the next row of data into the given row object.
	 *
	 * \param row A FlameRow instance to populate with values from the next row.
	 * \return True if there is another row to read after the present row.
	 */
	bool next(FlameRow& row);
};


class CSVReader : public Reader {
private:
	std::vector<std::vector<std::string>> m_data;
	std::string m_filename;
	int m_idx;
	bool m_transpose;
	int m_minWlCol;
	int m_maxWlCol;
	int m_headerRows;
	int m_idCol;

	void load();

	void transpose();

	void loadBandMap();

public:
	CSVReader(const std::string& filename, bool transpose, int headerRows, int minWlCol, int maxWlCol, int idCol);

	/**
	 * Attempt to guess the default properties for the document.
	 *
	 * The assumption is that the header will contain wavelengths, so there must be a contiguous
	 * sequence of floating point numbers. Each row will contain a sequence of floats, with
	 * zero or more columns of strings always at the same index. If a row contains a string field,
	 * the column header must not be a float.
	 *
	 * \param filename The CSV file.
	 * \param[out] transpose True if the document should be transposed.
	 * \param[out] header The number of header rows.
	 * \param[out] minCol The first data column.
	 * \param[out] maxCol The last data column.
	 * \param[out] idCol The likely identifier column, or -1 if there isn't one. Generally, the first non-numeric column.
	 */
	static void guessFileProperties(const std::string& filename, bool& transpose, int& header, int& minCol, int& maxCol, int& idCol);

	void reset();

	bool next(std::string& id, std::vector<double>& buf, int& cols, int& col, int& row);

};

} // reader
} // hlrg

#endif /* READER_HPP_ */
