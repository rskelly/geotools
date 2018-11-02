/*
 * raster.hpp
 *
 *  Created on: Sep 20, 2018
 *      Author: rob
 */

#ifndef INCLUDE_RASTER_HPP_
#define INCLUDE_RASTER_HPP_

#include <vector>
#include <string>

#include <gdal_priv.h>


namespace hlrg {


enum DataType {
	UInt16,
	Int16,
	UInt32,
	Int32,
	Float32,
	Float64,
	Byte
};

/**
 * Provides methods for reading rows out of a raster.
 * One might think this should be a Reader implementation, but its purpose
 * is just to read rows, not fit into spaces where a Reader would be appropriate.
 */
class Raster {
private:
	GDALDataset* m_ds;
	int m_bands;
	int m_cols;
	int m_rows;
	int m_row;

protected:

	GDALDataset* ds() const;

	GDALRasterBand* band(int b) const;

	char** metadata() const;

public:
	/**
	 * Construct a readable raster object using the given filename.
	 */
	Raster(const std::string& filename);

	/**
	 * Construct a writable raster object using the given filename.
	 */
	Raster(const std::string& filename, int cols, int rows, int bands, int srid, DataType type, Raster* parent = nullptr);

	/**
	 * Return the number of bands.
	 *
	 * @return The number of bands.
	 */
	int bands() const;

	/**
	 * Return the number of columns.
	 *
	 * @return The number of columns.
	 */
	int cols() const;


	int rows() const;


	/**
	 * Read all bands of the next row. This will write the pixel rows for each band
	 * end-to-end into the given buffer.
	 *
	 * @param buf A vector which contains the pixel row for each band, end to end.
	 * @return True if there is another row to be read. False otherwise.
	 */
	template <class T>
	bool next(std::vector<T>& buf);

	/**
	 * Write a row of data to the raster. This buffer contains the pixels of the row,
	 * with subsequent bands stored sequentially.
	 *
	 * @param buf The data to write.
	 * @param row The row index to write to.
	 * @return True if the write succeeds.
	 */
	template <class T>
	bool write(std::vector<T>& buf, int row);
	/**
	 * Get the row by zero-based index.
	 *
	 * @param buf A vector which contains the pixel row for each band, end to end.
	 * @param row The zero-based row index.
	 * @return True if the read was successful.
	 */
	template <class T>
	bool get(std::vector<T>& buf, int row);

	/**
	 * Reset the file pointer to the beginning.
	 */
	void reset();

	~Raster();
};


}


#endif /* INCLUDE_RASTER_HPP_ */
