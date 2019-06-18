/*
 * contrem_util.hpp
 *
 *  Created on: Jun 4, 2019
 *      Author: rob
 */

#ifndef INCLUDE_UTIL_HPP_
#define INCLUDE_UTIL_HPP_

#include <array>

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
 * Return the file descriptor of a temporary opened file.
 */
int tmpfile();

/**
 * Remove non-alphanumeric characters and replace with underscores.
 */
std::string sanitize(const std::string& str);

} // util
} // hlrg

#endif /* INCLUDE_UTIL_HPP_ */
