/*
 * contrem_util.cpp
 *
 *  Created on: Jun 4, 2019
 *      Author: rob
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>

#include <cstdio>
#include <algorithm>
#include <sstream>
#include <regex>

#include <gdal_priv.h>

#include "util.hpp"
#include "reader.hpp"

using namespace hlrg::util;

FileType hlrg::util::getFileType(const std::string& filename) {
	std::string ext;
	{
		size_t p = filename.find('.');
		if(p < std::string::npos) {
			std::string ext0 = filename.substr(p);
			std::transform(ext0.begin(), ext0.end(), std::back_inserter(ext), ::tolower);
		}
	}
	if(ext == ".csv" || ext == ".txt") {
		return FileType::CSV;
	} else if(ext == ".roi") {
		return FileType::ROI;
	} else {
		GDALAllRegister();
		GDALDataset* ds = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_READONLY, 0, 0, 0));
		if(ds) {
			std::string drv(ds->GetDriverName());
			FileType type = FileType::Unknown;
			if(drv == "GTiff") {
				type = FileType::GTiff;
			} else if(drv == "ENVI") {
				type = FileType::ENVI;
			} else if(drv == "ESRI Shapefile") {
				type = FileType::SHP;
			} else if(drv == "SQLite") {
				type = FileType::SQLITE;
			}
			GDALClose(ds);
			return type;
		}
	}
	return FileType::Unknown;
}

std::string hlrg::util::fileTypeAsString(FileType type) {
	switch(type) {
	case FileType::GTiff: return "GTiff";
	case FileType::ENVI: return "ENVI";
	case FileType::ROI: return "ENVI ROI";
	case FileType::SHP: return "Shapefile";
	case FileType::CSV: return "CSV";
	default: return "";
	}
}

FileType hlrg::util::fileTypeFromString(const std::string& type) {
	if(type == "GTiff") {
		return FileType::GTiff;
	} else if(type == "ENVI") {
		return FileType::ENVI;
	} else if(type == "ENVI ROI") {
		return FileType::ROI;
	} else if(type == "Shapefile" || type == "ESRI Shapefile") {
		return FileType::SHP;
	} else if(type == "CSV") {
		return FileType::CSV;
	} else {
		return FileType::Unknown;
	}
}

std::string hlrg::util::normMethodAsString(NormMethod method) {
	switch(method) {
	case NormMethod::ConvexHull:
		return "Convex Hull";
	case NormMethod::ConvexHullLongestSeg:
		return "Convex Hull, Longest Segment";
	case NormMethod::Line:
		return "Line";
	case NormMethod::Unknown:
	default:
		return "Unknown";
	}
}

NormMethod hlrg::util::normMethodFromString(const std::string& method) {
	if(method == normMethodAsString(NormMethod::ConvexHull)) {
		return NormMethod::ConvexHull;
	} else if(method == normMethodAsString(NormMethod::ConvexHullLongestSeg)) {
		return NormMethod::ConvexHullLongestSeg;
	} else if(method == normMethodAsString(NormMethod::Line)) {
		return NormMethod::Line;
	} else {
		return NormMethod::Unknown;
	}
}

bool hlrg::util::isnonzero(const double& v) {
	return v != 0;
}

bool hlrg::util::isdir(const std::string& path) {
	struct stat st;
	if(!stat(path.c_str(), &st))
		return S_ISDIR(st.st_mode);
	return false;
}

bool hlrg::util::isfile(const std::string& path) {
	struct stat st;
	if(!stat(path.c_str(), &st))
		return S_ISREG(st.st_mode);
	return false;
}

bool hlrg::util::rem(const std::string& dir) {
	return !unlink(dir.c_str());
}

bool hlrg::util::makedir(const std::string& filename) {
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

int hlrg::util::tmpfile() {
	char tpl[] = {"geo_util_XXXXXX"};
	return mkstemp(tpl);
}

std::string hlrg::util::sanitize(const std::string& str) {
	std::regex repl("([^0-9A-Za-z]+)");
	std::stringstream ss;
	std::regex_replace(std::ostreambuf_iterator<char>(ss), str.begin(), str.end(), repl, "_");
	return ss.str();
}
