/*
 * contrem_util.cpp
 *
 *  Created on: Jun 4, 2019
 *      Author: rob
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <sstream>
#include <regex>

#include <gdal_priv.h>

#include "reader.hpp"
#include "contrem_util.hpp"

using namespace hlrg;

FileType hlrg::getFileType(const std::string& filename) {
	std::string ext;
	{
		size_t p = filename.find('.');
		if(p < std::string::npos) {
			std::string ext0 = filename.substr(p);
			std::transform(ext0.begin(), ext0.end(), std::back_inserter(ext), ::tolower);
		}
	}
	if(ext == ".csv") {
		return CSV;
	} else if(ext == ".roi") {
		return ROI;
	} else {
		GDALAllRegister();
		GDALDataset* ds = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_READONLY, 0, 0, 0));
		if(ds) {
			std::string drv(ds->GetDriverName());
			FileType type = UnknownFileType;
			if(drv == "GTiff") {
				type = GTiff;
			} else if(drv == "ENVI") {
				type = ENVI;
			} else if(drv == "ESRI Shapefile") {
				type = SHP;
			} else if(drv == "SQLite") {
				type = SQLITE;
			}
			GDALClose(ds);
			return type;
		}
	}
	return UnknownFileType;
}

/**
 * Return a map containing pairs where the int is the 1-based band index,
 * and the float is the wavelength. Attempts to load from raster metadata
 * or table header. If these fail, will attempt to load from first column
 * of presumably transposed table.
 */
std::map<int, double> hlrg::loadWavelengths(const Contrem& contrem) {
	std::map<int, double> map;
	switch(getFileType(contrem.spectra)) {
	case GTiff:
	case ENVI:
		{
			GDALReader rdr(contrem.spectra);
			for(const auto& it : rdr.getBandMap())
				map[it.second] = (double) it.first / WL_SCALE;
		}
		break;
	case CSV:
		{
			CSVReader rdr(contrem.spectra, contrem.wlTranspose, contrem.wlHeaderRows, contrem.wlMinCol, contrem.wlMaxCol, contrem.wlIDCol);
			for(const auto& it : rdr.getBandMap())
				map[it.second] = (double) it.first / WL_SCALE;
		}
		break;
	case SHP:
	case ROI:
	default:
		throw std::runtime_error("Invalid file type: " + contrem.spectra);
	}
	return map;
}

std::string hlrg::fileTypeAsString(FileType type) {
	switch(type) {
	case GTiff: return "GTiff";
	case ENVI: return "ENVI";
	case ROI: return "ENVI ROI";
	case SHP: return "Shapefile";
	case CSV: return "CSV";
	default: return "";
	}
}

hlrg::FileType hlrg::fileTypeFromString(const std::string& type) {
	if(type == "GTiff") {
		return GTiff;
	} else if(type == "ENVI") {
		return ENVI;
	} else if(type == "ENVI ROI") {
		return ROI;
	} else if(type == "Shapefile" || type == "ESRI Shapefile") {
		return SHP;
	} else if(type == "CSV") {
		return CSV;
	} else {
		return UnknownFileType;
	}
}

std::string hlrg::normMethodAsString(hlrg::NormMethod method) {
	switch(method) {
	case ConvexHull:
		return "Convex Hull";
	case ConvexHullLongestSeg:
		return "Convex Hull, Longest Segment";
	case Line:
		return "Line";
	case UnknownNormMethod:
	default:
		return "Unknown";
	}
}

hlrg::NormMethod hlrg::normMethodFromString(const std::string& method) {
	if(method == normMethodAsString(ConvexHull)) {
		return ConvexHull;
	} else if(method == normMethodAsString(ConvexHullLongestSeg)) {
		return ConvexHullLongestSeg;
	} else if(method == normMethodAsString(Line)) {
		return Line;
	} else {
		return UnknownNormMethod;
	}
}

bool hlrg::isnonzero(const double& v) {
	return v != 0;
}

bool hlrg::isdir(const std::string& path) {
	struct stat st;
	if(!stat(path.c_str(), &st))
		return S_ISDIR(st.st_mode);
	return false;
}

bool hlrg::isfile(const std::string& path) {
	struct stat st;
	if(!stat(path.c_str(), &st))
		return S_ISREG(st.st_mode);
	return false;
}

bool hlrg::rem(const std::string& dir) {
	return !unlink(dir.c_str());
}

bool hlrg::makedir(const std::string& filename) {
	std::stringstream path(filename);
	std::stringstream inter;
	std::string part, current;
	if(filename[0] == '/')
		inter << '/';
	while(std::getline(path, part, '/')) {
		inter << part;
		current = inter.str();
		if(!isdir(current) && !isfile(current)) {
			if(mkdir(current.c_str(), 0755))
				return false;
		}
		if(!part.empty())
			inter << '/';
	}
	return true;
}

std::string hlrg::sanitize(const std::string& str) {
	std::regex repl("([^0-9A-Za-z]+)");
	std::stringstream ss;
	std::regex_replace(std::ostreambuf_iterator<char>(ss), str.begin(), str.end(), repl, "_");
	return ss.str();
}
