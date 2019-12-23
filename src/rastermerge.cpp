/*
 * rastermerge.cpp
 *
 *  Created on: Dec 19, 2019
 *      Author: rob
 */

#include <string>
#include <vector>
#include <iostream>
#include <list>
#include <thread>
#include <mutex>

#include <gdal_priv.h>

class Ctx {
public:
	std::vector<float> data;
	int cols;
	int rows;
	double trans[6];
	float nd;
};

bool loadRaster(const std::string& file, int band, Ctx& data) {

	GDALAllRegister();

	GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(file.c_str(), GA_ReadOnly));
	if(!ds)
		return false;
	ds->GetGeoTransform(data.trans);
	data.cols = ds->GetRasterXSize();
	data.rows = ds->GetRasterYSize();
	GDALRasterBand* bnd = ds->GetRasterBand(band);
	data.nd = bnd->GetNoDataValue();
	data.data.resize(data.cols * data.rows);
	if(CE_None != bnd->RasterIO(GF_Read, 0, 0, data.cols, data.rows, data.data.data(), data.cols, data.rows, GDT_Float32, 0, 0, 0)) {
		GDALClose(ds);
		return false;
	}
	GDALClose(ds);
	return true;
}

bool saveRaster(const std::string& tpl, const std::string& file, std::vector<float>& data) {
	GDALDataset* dst = static_cast<GDALDataset*>(GDALOpen(tpl.c_str(), GA_ReadOnly));
	if(!dst)
		return false;
	int cols = dst->GetRasterXSize();
	int rows = dst->GetRasterYSize();
	double trans[6];
	dst->GetGeoTransform(trans);
	GDALDriver* drv = dst->GetDriver();
	GDALDataset* ds = static_cast<GDALDataset*>(drv->Create(file.c_str(), cols, rows, 1, GDT_Float32, 0));
	if(!ds) {
		GDALClose(dst);
		return false;
	}
	ds->SetProjection(dst->GetProjectionRef());
	ds->SetGeoTransform(trans);
	ds->GetRasterBand(1)->SetNoDataValue(dst->GetRasterBand(1)->GetNoDataValue());
	GDALClose(dst);
	if(CE_None != ds->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, cols, rows, data.data(), cols, rows, GDT_Float32, 0, 0, 0)) {
		GDALClose(ds);
		return false;
	}
	GDALClose(ds);
	return true;
}

int toCol(double x, double* trans) {
	return (int) (x - trans[0]) / trans[1];
}

int toRow(double y, double* trans) {
	return (int) (y - trans[3]) / trans[5];
}

double toX(int col, double* trans) {
	return (col * trans[1]) + trans[0] + trans[1] * 0.5;
}

double toY(int row, double* trans) {
	return (row * trans[5]) + trans[3] + trans[5] * 0.5;
}



bool smooth(std::vector<bool>& filled, std::vector<float>& src, std::vector<float>& dst, int col, int row, int cols, int rows) {
	if(filled[row * cols + col])
		return false;
	float t = 0, w = 0;
	int n = 0;
	double v;
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {
			if(!std::isnan((v = src[r * cols + c]))) {
				float w0 = c == col && r == row ? 1.0 : 1.0 / (std::pow(c - col, 2) + std::pow(r - row, 2));
				t += v * w0;
				w += w0;
				++n;
			}
		}
	}
	if(n && w > 0) {
		dst[row * cols + col] = t / w;
		filled[row * cols + col] = true;
		return true;
	}
	return false;
}

void processIDW(std::list<int>* rowq, std::mutex* qmtx, Ctx* src, Ctx* dst, std::mutex* dmtx, int size) {
	int row;
	while(!rowq->empty()) {
		{
			std::lock_guard<std::mutex> lk(*qmtx);
			if(!rowq->empty()) {
				row = rowq->front();
				rowq->pop_front();
				std::cout << "Row " << row << " of " << src->rows << "\n";
			} else {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
				continue;
			}
		}
		double v0;
		for(int col = 0; col < src->cols; ++col) {
			float s = 0;
			float w = 0;
			bool halt = false;
			for(int r = -size / 2; !halt && r < size / 2 + 1; ++r) {
				for(int c = -size / 2; c < size / 2 + 1; ++c) {
					int cc = col + c;
					int rr = row + r;
					if(cc < 0 || rr < 0 || cc >= src->cols || rr >= src->rows)
						continue;
					if(!std::isnan((v0 = src->data[rr * src->cols + cc]))) {
						float d = (float) (c * c + r * r);
						if(d == 0) {
							s = v0;
							w = 1;
							break;
						} else {
							float w0 = 1.0 / d;
							s += v0 * w0;
							w += w0;
						}
					}
				}
			}
			if(w > 0) {
				std::lock_guard<std::mutex> lk(*dmtx);
				dst->data[row * src->cols + col] = s / w;
			}
		}
	}
}

void processCos(std::list<int>* rowq, std::mutex* qmtx, Ctx* src, Ctx* dst, std::mutex* dmtx, int size, float* cos) {
	float rad2 = std::pow(size / 2.0, 2.0) / 1000;	// Because the cosine lookup has 1000 elements.
	int row;
	while(!rowq->empty()) {
		{
			std::lock_guard<std::mutex> lk(*qmtx);
			if(!rowq->empty()) {
				row = rowq->front();
				rowq->pop_front();
				std::cout << "Row " << row << " of " << src->rows << "\n";
			} else {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
				continue;
			}
		}

		float v0;
		for(int col = 0; col < src->cols; ++col) {
			float s = 0;
			float w = 0;
			bool halt = false;
			for(int r = -size / 2; !halt && r < size / 2 + 1; ++r) {
				for(int c = -size / 2; c < size / 2 + 1; ++c) {
					int cc = col + c;
					int rr = row + r;
					if(cc < 0 || rr < 0 || cc >= src->cols || rr >= src->rows)
						continue;
					if(!std::isnan((v0 = src->data[rr * src->cols + cc]))) {
						int d = (int) std::min(1000.0f, (float) (c * c + r * r) / rad2);
						//std::cout << d << " " << c << " " << r << " " << (c * c + r * r) << " " << rad2 << "\n";
						if(d < 1000){
							float w0 = cos[d];
							s += v0 * w0;
							w += 1;
						}
					}
				}
			}
			if(w > 0) {
				std::lock_guard<std::mutex> lk(*dmtx);
				dst->data[row * src->cols + col] = s / w;
			}
		}
	}
}


void processGauss(std::list<int>* rowq, std::mutex* qmtx, Ctx* src, Ctx* dst, std::mutex* dmtx, int size, float sigma) {
	int row;
	while(!rowq->empty()) {
		{
			std::lock_guard<std::mutex> lk(*qmtx);
			if(!rowq->empty()) {
				row = rowq->front();
				rowq->pop_front();
				std::cout << "Row " << row << " of " << src->rows << "\n";
			} else {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
				continue;
			}
		}

		float gain = 0.5;
		float v0;
		for(int col = 0; col < src->cols; ++col) {
			float s = 0;
			float w = 0;
			bool halt = false;
			for(int r = -size / 2; !halt && r < size / 2 + 1; ++r) {
				for(int c = -size / 2; c < size / 2 + 1; ++c) {
					int cc = col + c;
					int rr = row + r;
					if(cc < 0 || rr < 0 || cc >= src->cols || rr >= src->rows)
						continue;
					if(!std::isnan((v0 = src->data[rr * src->cols + cc]))) {
						float w0 = std::exp(-0.5 * (c * c + r * r) / (sigma * sigma));
						//std::cout << d << " " << c << " " << r << " " << (c * c + r * r) << " " << rad2 << "\n";
						s += v0 * w0;
						w += gain;
					}
				}
			}
			if(w > 0) {
				std::lock_guard<std::mutex> lk(*dmtx);
				dst->data[row * src->cols + col] = s / w;
			}
		}
	}
}

int main(int argc, char** argv) {

	if(argc < 6) {
		std::cerr << "Usage: rastermerge [options] <input file 1> <input band 1> <input file 2> <input band 2> <output file>\n"
				<< " -s <size>    The size of the window in pixels.\n"
				<< " -t <threads> The number of threads.\n"
				<< " -d           Preserve the difference raster in /tmp/diff.tif\n";
		return 1;
	}

	std::vector<std::string> files;
	std::vector<int> bands;
	int size = 501;
	int tcount = 4;
	bool diff = false;

	for(int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-s") {
			size = atoi(argv[++i]);
			if(size % 2 == 0)
				++size;
		} else if(arg == "-d") {
			diff = true;
		} else if(arg == "-t") {
			tcount = atoi(argv[++i]);
			if(tcount < 1)
				tcount = 1;
		} else {
			files.push_back(arg);
			if(i < argc - 1)
				bands.push_back(atoi(argv[++i]));
		}
	}

	std::string file1 = files[0];
	int band1 = bands[0];
	std::string file2 = files[1];
	int band2 = bands[1];

	std::string outfile = files[2];

	std::cout << "Matching: " << file2 << " (" << band2 << ") to " << file1 << " (" << band1 << ") --> " << outfile << "; size: " << size << "\n";

	Ctx data2;
	loadRaster(file2, band2, data2);

	Ctx src = data2;
	std::vector<bool> filled(src.cols * src.rows);
	std::fill(src.data.begin(), src.data.end(), std::nan(""));
	std::fill(filled.begin(), filled.end(), false);

	{
		Ctx data1;
		loadRaster(file1, band1, data1);

		for(int row2 = 0; row2 < data2.rows; ++row2) {
			for(int col2 = 0; col2 < data2.cols; ++col2) {
				double v2 = data2.data[row2 * data2.cols + col2];
				if(v2 != data2.nd) {
					double x = toX(col2, data2.trans);
					double y = toY(row2, data2.trans);
					int col1 = toCol(x, data1.trans);
					int row1 = toRow(y, data1.trans);
					if(!(col1 < 0 || row1 < 0 || col1 >= data1.cols || row1 >= data1.rows)) {
						double v1 = data1.data[row1 * data1.cols + col1];
						if(v1 != data1.nd)
							src.data[row2 * data2.cols + col2] = v1 - v2;
					}
				}
			}
		}
	}

	if(diff)
		saveRaster(file2, "/tmp/diff.tif", src.data);

	Ctx dst = src;
	std::fill(dst.data.begin(), dst.data.end(), 0);

	{

		float cos[1001];
		for(int i = 0; i <= 1000; ++i)
			cos[i] = std::cos((float) i / 1000 * M_PI) / 2.0 + 0.5;

		std::list<int> rowq;
		for(int row = 0; row < data2.rows; ++row)
			rowq.push_back(row);

		std::mutex qmtx;
		std::mutex dmtx;
		std::vector<std::thread> threads;
		for(int i = 0; i < tcount; ++i) {
			threads.emplace_back(processGauss, &rowq, &qmtx, &src, &dst, &dmtx, size, (float) size * 0.25);
		}
		for(int i = 0; i < tcount; ++i) {
			if(threads[i].joinable())
				threads[i].join();
		}
	}

	saveRaster(file2, outfile, dst.data);

	return 0;
}

