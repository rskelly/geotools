/*
 * contrem_util.hpp
 *
 *  Created on: Jun 4, 2019
 *      Author: rob
 */

#ifndef INCLUDE_CONTREM_UTIL_HPP_
#define INCLUDE_CONTREM_UTIL_HPP_

#include <array>

#include "contrem.hpp"

namespace hlrg {

	constexpr std::array<FileType, 3> OUTPUT_TYPES = {GTiff, ENVI, CSV};			///<! Allowed output types for results.

	FileType getFileType(const std::string& filename);

	/**
	 * Return a map containing pairs where the int is the 1-based band index,
	 * and the float is the wavelength. Attempts to load from raster metadata
	 * or table header. If these fail, will attempt to load from first column
	 * of presumably transposed table.
	 */
	std::map<int, double> loadWavelengths(const Contrem& contrem);

	std::string fileTypeAsString(hlrg::FileType type);

	hlrg::FileType fileTypeFromString(const std::string& type);

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


} // hlrg

#endif /* INCLUDE_CONTREM_UTIL_HPP_ */
