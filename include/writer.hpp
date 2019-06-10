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

namespace hlrg {

/**
 * Enumeration containing common data types.
 */
enum DataType {
	Byte,
	Int32,
	Float32
};


/**
 * An abstract class used to write values and statistics in a standard way.
 */
class Writer {
public:
	/**
	 * Write the given buffer using the grid coordinates.
	 *
	 * \param buf A vector containing values to write.
	 * \param col The column offset to write to.
	 * \param row The row offset to write to.
	 * \param cols The number of columns to write.
	 * \param rows The number of rows to write.
	 * \param buffSizeX The size of the buffer in the x dimension.
	 * \param buffSizeY The size of the buffer in the y dimension.
	 */
	virtual bool write(const std::vector<double>& buf, int col, int row, int cols, int rows, int bufSizeX = 0, int bufSizeY = 0) = 0;

	/**
	 * Write the given buffer using the grid coordinates.
	 *
	 * \param buf A vector containing values to write.
	 * \param col The column offset to write to.
	 * \param row The row offset to write to.
	 * \param cols The number of columns to write.
	 * \param rows The number of rows to write.
	 * \param buffSizeX The size of the buffer in the x dimension.
	 * \param buffSizeY The size of the buffer in the y dimension.
	 */
	virtual bool write(const std::vector<int>& buf, int col, int row, int cols, int rows, int bufSizeX = 0, int bufSizeY = 0) = 0;

	/**
	 * Compute stats and write them to the file with the given filename.
	 * The names list is optional, and used for naming the individual statistics.
	 *
	 * \param filename The output filename.
	 * \param names A list of names for the statistics.
	 * \return True if successful.
	 */
	virtual bool writeStats(const std::string& filename, const std::vector<std::string>& names = {}) = 0;

	/**
	 * Fill the dataset with the given value.
	 *
	 * \param v The value to fill with.
	 */
	virtual void fill(double v) = 0;

	virtual ~Writer() {}
};

/**
 * An implementation of Writer that can write to GDAL datasources.
 */
class GDALWriter : public Writer {
private:
	GDALDataset* m_ds;
	int m_bands;
	int m_cols;
	int m_rows;

public:

	/**
	 * Construct a GDALWriter.
	 *
	 * \param filename The output file name.
	 * \param driver The file format driver name.
	 * \param cols The size of the output in columns.
	 * \param rows The size of the output in rows.
	 * \param bands The number of bands.
	 * \param wavelengths A list of the wavelengths corresponding to each band.
	 * \param bandNames A list of the names of the bands corresponding to each band.
	 * \param type The data type of the output.
	 * \param interleave The interleaving format.
	 * \param unit The wavelength units.
	 */
	GDALWriter(const std::string& filename, const std::string& driver, int cols, int rows, int bands,
			const std::vector<double>& wavelengths = {}, const std::vector<std::string>& bandNames = {}, char** meta = nullptr,
			DataType type = DataType::Float32, const std::string& interleave = "BIL", const std::string& unit = "nm");

	bool write(const std::vector<double>& buf, int col, int row, int cols, int rows, int bufSizeX = 0, int bufSizeY = 0);

	bool write(const std::vector<int>& buf, int col, int row, int cols, int rows, int bufSizeX = 0, int bufSizeY = 0);

	bool writeStats(const std::string& filename, const std::vector<std::string>& names = {});

	/**
	 * Fill the raster with the given value.
	 *
	 * \param v The value to fill with.
	 */
	void fill(double v);

	/**
	 * Fill the raster with the given value.
	 *
	 * \param v The value to fill with.
	 */
	void fill(int v);

	~GDALWriter();
};

} // hlrg

#endif /* WRITER_HPP_ */
