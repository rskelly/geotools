/*
 * contrem_util.hpp
 *
 *  Created on: Jun 4, 2019
 *      Author: rob
 */

#ifndef INCLUDE_UTIL_HPP_
#define INCLUDE_UTIL_HPP_

#include <array>
#include <cstring>

#include <gdal_priv.h>

namespace hlrg {
namespace util {

/**
 * Input and output file types.
 */
enum class FileType {
	GTiff,
	ENVI,
	ROI,
	SHP,
	CSV,
	SQLITE,
	Unknown
};

/**
 * Enumeration containing common data types.
 */
enum class DataType {
	Byte,
	Int32,
	Float32
};

/**
 * Normalization methods.
 */
enum class NormMethod {
	ConvexHull,
	ConvexHullLongestSeg,
	Line,
	Unknown
};

constexpr std::array<FileType, 3> OUTPUT_TYPES = {FileType::GTiff, FileType::ENVI, FileType::CSV};			///<! Allowed output types for results.

constexpr std::array<NormMethod, 3> NORM_METHODS = {NormMethod::ConvexHull, NormMethod::ConvexHullLongestSeg, NormMethod::Line};

/**
 * A class that when instantiated creates a temporary file
 * whose lifecycle is automatically managed. When the class
 * destructs, the file is deleted. Maintains the filename
 * and file descriptor.
 */
class TmpFile {
public:
	std::string filename;	///<! The filename of the temporary file.
	int fd;					///<! The file descriptor.
	size_t size;			///<! The file size.

	/**
	 * Create a file with the given size. The contents of the file are not defined.
	 *
	 * \param size The file size.
	 */
	TmpFile(size_t size = 0);

	/**
	 * Resize the file.
	 *
	 * \param size The size.
	 */
	void resize(size_t size);

	/**
	 * Close the file.
	 */
	void close();

	/**
	 * Destroy. Closes and deletes the file.
	 */
	~TmpFile();
};


FileType getFileType(const std::string& filename);

std::string fileTypeAsString(FileType type);

FileType fileTypeFromString(const std::string& type);

NormMethod normMethodFromString(const std::string& method);

std::string normMethodAsString(NormMethod method);

bool isnonzero(const double& v);

/**
 * Return true if it's a dir and it exists.
 */
bool isdir(const std::string& path);

/**
 * Return true if it's a file and it exists.
 */
bool isfile(const std::string& path);

/**
 * Remove the directory or file.
 */
bool rem(const std::string& dir);

/**
 * Recursively make the directory.
 */
bool makedir(const std::string& filename);

/**
 * Remove non-alphanumeric characters and replace with underscores.
 */
std::string sanitize(const std::string& str);

/**
 * Return the byte size of the given GDAL type.
 *
 * \param The GDAL type.
 * \return The type size.
 */
int gdalTypeSize(GDALDataType type);

/**
 * Convert the given char buffer to a buffer of the templated type.
 *
 * \param[in] The source buffer.
 * \param[out] The output buffer.
 */
template <class T>
void inline convertBuffer(std::vector<char>& raw, std::vector<T>& buf) {
	buf.reserve(raw.size() / sizeof(T));
	std::memcpy(buf.data(), raw.data(), raw.size());
}

/**
 * Convert the given raw char buffer to the typed output buffer
 * using the GDALDataType as a guide.
 *
 * \param type the GDALDataType.
 * \param[in] rawBuf The raw character buffer.
 * \param[out] buf The output buffer.
 */
template <class T>
void convertBuffer(GDALDataType type, std::vector<char>& rawBuf, std::vector<T>& buf) {

	switch(type) {
	case GDT_Float32:
	{
		std::vector<float> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_Int32:
	{
		std::vector<int32_t> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_UInt32:
	{
		std::vector<uint32_t> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_Int16:
	{
		std::vector<int16_t> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_UInt16:
	{
		std::vector<uint16_t> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_Float64:
	{
		std::vector<double> tmp;
		convertBuffer(rawBuf, tmp);
		std::copy(tmp.begin(), tmp.end(), buf.begin());
		break;
	}
	case GDT_Byte:
		std::copy(rawBuf.begin(), rawBuf.end(), buf.begin());
		break;
	default:
		throw std::runtime_error("Unknown GDAL data type: " + std::to_string((int) type));
	}
}


} // util
} // hlrg

#endif /* INCLUDE_UTIL_HPP_ */
