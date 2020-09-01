/*
 * map2raw.cpp
 *
 * Converts a GIS raster to a raw format, in 8 or 16 bits with no header.
 * Used for height maps in Unity, etc.
 *
 *  Created on: Aug 23, 2020
 *      Author: rob
 */

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include <gdal_priv.h>

int main(int, char** argv) {

	int p = 1;
	std::string infile(argv[p++]);
	int band = atoi(argv[p++]);
	int bits = atoi(argv[p++]);
	//float minx = atof(argv[p++]);
	//float miny = atof(argv[p++]);
	//float maxx = atof(argv[p++]);
	//float maxy = atof(argv[p++]);
	std::string outfile(argv[p++]);

	GDALAllRegister();
	GDALDataset* ds = (GDALDataset*) GDALOpen(infile.c_str(), GA_ReadOnly);
	double trans[6];
	ds->GetGeoTransform(trans);
	int cols = ds->GetRasterXSize();
	int rows = ds->GetRasterYSize();
	GDALRasterBand* bnd = ds->GetRasterBand(band);
	//GDALDataType type = bnd->GetRasterDataType();
	double nd = bnd->GetNoDataValue();

	std::vector<float> buf(cols * rows);
	if(CE_None != bnd->RasterIO(GF_Read, 0, 0, cols, rows, buf.data(), cols, rows, GDT_Float32, 0, 0, 0)) {
		std::cerr << "Failed to read raster\n";
		return 1;
	}

	float min = std::numeric_limits<float>::max();
	float max = std::numeric_limits<float>::lowest();
	for(size_t i = 0; i < buf.size(); ++i) {
		if(buf[i] < min) min = buf[i];
		if(buf[i] > max) max = buf[i];
	}

	int extra_cols = 0;
	int extra_rows = 0;
	int dif;
	if((dif = (cols + extra_cols) % 4096) > 0) {
		extra_cols += dif;
	}
	if((dif = (rows + extra_rows) % 4096) > 0) {
		extra_rows += dif;
	}
	int bcols = (cols + extra_cols) / 4096 + 1;
	int brows = (rows + extra_rows) / 4096 + 1;
	for(int bcol = 0; bcol < bcols; ++bcol) {
		for(int brow = 0; brow < brows; ++brow) {
			std::string base = outfile.substr(0, outfile.find("."));
			std::stringstream outf;
			outf << base << "_" << bcol << "_" << brow << outfile.substr(outfile.find("."));
			std::string outff = outf.str();
			std::ofstream out(outff, std::ios::binary);

			for(int col = bcol * 4096; col < (bcol + 1) * 4096; ++col) {
				for(int row = brow * 4096; row < (brow + 1) * 4096; ++row) {
					float v = col >= cols || row >= rows ? std::nan("") : buf[row * cols + col];
					if((double) v == nd)
						v = std::nan("");
					if(bits == 8) {
						uint8_t b = std::isnan(v) ? 0 : (uint8_t) (v - min) / (max - min) * 255;
						out << b;
					} else {
						uint16_t b = std::isnan(v) ? 0 : (uint8_t) (v - min) / (max - min) * 65535;
						out << b;
					}
				}
			}
		}
	}
}


