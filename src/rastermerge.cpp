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

#include "ds/mqtree.hpp"

class Ctx {
public:
	std::vector<float> data;
	int cols;
	int rows;
	double trans[6];
	float nd;
	std::string projection;

	Ctx& operator=(const Ctx& other) {
		cols = other.cols;
		rows = other.rows;
		nd = other.nd;
		projection = other.projection;
		data.resize(other.data.size());
		for(int i = 0; i < 6; ++i)
			trans[i] = other.trans[i];
		return *this;
	}
};

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

bool loadRaster(const std::string& file, int band, Ctx& data) {

	GDALAllRegister();

	GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(file.c_str(), GA_ReadOnly));
	if(!ds)
		return false;
	ds->GetGeoTransform(data.trans);
	data.cols = ds->GetRasterXSize();
	data.rows = ds->GetRasterYSize();
	ds->GetGeoTransform(data.trans);
	GDALRasterBand* bnd = ds->GetRasterBand(band);
	data.nd = bnd->GetNoDataValue();
	data.data.resize(data.cols * data.rows);
	data.projection = ds->GetProjectionRef();
	if(CE_None != bnd->RasterIO(GF_Read, 0, 0, data.cols, data.rows, data.data.data(), data.cols, data.rows, GDT_Float32, 0, 0, 0)) {
		GDALClose(ds);
		return false;
	}
	GDALClose(ds);
	return true;
}

bool loadRasters(const std::vector<std::string>& files, std::vector<int> bands, int count, int resample, Ctx& data) {

	GDALAllRegister();

	double trans[6];
	double minx = 99999999.0, miny = 999999999.0;
	double maxx = -99999999.0, maxy = -99999999.0;
	int cols, rows;

	for(size_t i = 0; i < (size_t) count; ++i) {
		std::string file = files[i];
		int band = bands[i];

		GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(file.c_str(), GA_ReadOnly));
		if(!ds)
			return false;
		ds->GetGeoTransform(trans);
		cols = ds->GetRasterXSize();
		rows = ds->GetRasterYSize();
		data.nd = ds->GetRasterBand(band)->GetNoDataValue();
		data.projection = ds->GetProjectionRef();

		double bounds[4] = {
			trans[1] > 0 ? trans[0] : trans[0] + cols * trans[1],
			trans[5] > 0 ? trans[3] : trans[3] + rows * trans[5],
			trans[1] < 0 ? trans[0] : trans[0] + cols * trans[1],
			trans[5] < 0 ? trans[3] : trans[3] + rows * trans[5],
		};

		if(bounds[0] < minx) minx = bounds[0];
		if(bounds[1] < miny) miny = bounds[1];
		if(bounds[2] > maxy) maxx = bounds[2];
		if(bounds[3] > maxy) maxy = bounds[3];

		GDALClose(ds);
	}

	trans[0] = trans[1] > 0 ? minx : maxx;
	trans[3] = trans[5] > 0 ? miny : maxy;
	trans[1] *= resample;
	trans[5] *= resample;
	data.cols = (int) std::ceil((maxx - minx) / std::abs(trans[1]));
	data.rows = (int) std::ceil((maxy - miny) / std::abs(trans[5]));
	data.data.resize(data.cols * data.rows);
	std::fill(data.data.begin(), data.data.end(), std::nan(""));
	for(int i = 0; i < 6; ++i)
		data.trans[i] = trans[i];

	std::vector<float> buf;

	for(size_t i = 0; i < (size_t) count; ++i) {
		std::string file = files[i];
		int band = bands[i];

		GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(file.c_str(), GA_ReadOnly));
		if(!ds)
			return false;
		ds->GetGeoTransform(trans);
		cols = ds->GetRasterXSize();
		rows = ds->GetRasterYSize();
		ds->GetGeoTransform(trans);
		buf.resize((cols / resample) * (rows / resample));
		GDALRasterBand* bnd = ds->GetRasterBand(band);
		if(CE_None != bnd->RasterIO(GF_Read, 0, 0, cols, rows, buf.data(), cols / resample, rows / resample, GDT_Float32, 0, 0, 0)) {
			GDALClose(ds);
			return false;
		}
		trans[1] *= resample;
		trans[5] *= resample;
		float nd = bnd->GetNoDataValue();
		for(int r = 0; r < rows / resample; ++r) {
			for(int c = 0; c < cols / resample; ++c) {
				float x = toX(c, trans);
				float y = toY(r, trans);
				int cc = toCol(x, data.trans);
				int rr = toRow(y, data.trans);
				if(!(cc < 0 || rr < 0 || cc >= data.cols || rr >= data.rows)) {
					double v = buf[r * (cols / resample) + c];
					if(v != nd)
						data.data[rr * data.cols + cc] = v;
				}
			}
		}

		GDALClose(ds);
	}
	return true;
}

bool saveRaster(const std::string& file, Ctx& data) {
	GDALDriverManager* dm = GetGDALDriverManager();
	GDALDriver* drv = dm->GetDriverByName("GTiff");
	GDALDataset* ds = drv->Create(file.c_str(), data.cols, data.rows, 1, GDT_Float32, 0);
	if(!ds)
		return false;
	ds->SetProjection(data.projection.c_str());
	ds->SetGeoTransform(data.trans);
	ds->GetRasterBand(1)->SetNoDataValue(data.nd);
	if(CE_None != ds->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, data.cols, data.rows, data.data.data(), data.cols, data.rows, GDT_Float32, 0, 0, 0)) {
		GDALClose(ds);
		return false;
	}
	GDALClose(ds);
	return true;
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

void processDW(std::list<int>* rowq, std::mutex* qmtx, Ctx* src, Ctx* dst, std::mutex* dmtx, int size) {
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
						float w0 = std::max(0.0, std::min(1.0, d / std::pow((float) size / 2.0, 2.0)));
						s += v0 * w0;
						w += w0;
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
		std::cerr << "Usage: rastermerge [options] <<anchor file 1> <anchor band 1> [<anchor file 2> <anchor band 2> [...]]> <target file 2> <target band 2> <output file>\n"
				<< " -s <size>          The size of the window in pixels.\n"
				<< " -t <threads>       The number of threads.\n"
				<< " -d                 Preserve the difference raster in /tmp/diff.tif\n"
				<< " -m <method>        The method: idw, dw, gauss, cosine. Default IDW.\n"
				<< " -k <mask> <band>   A mask file. Pixel value 1 is kept.\n"
				<< " -r <resample>      Resample the anchor rasters. 1 is native resolution; 2 is half; 4 is quarter, etc.\n";
		return 1;
	}

	std::vector<std::string> files;
	std::vector<int> bands;
	float radius = 100;
	int tcount = 4;
	bool diff = false;
	std::string method = "idw";
	int resample = 1;
	std::string outfile;
	std::string maskfile;
	int maskband = 0;

	for(int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-s") {
			radius = atof(argv[++i]);
		} else if(arg == "-d") {
			diff = true;
		} else if(arg == "-r") {
			resample = atoi(argv[++i]);
		} else if(arg == "-k") {
			maskfile = argv[++i];
			maskband = atoi(argv[++i]);
		} else if(arg == "-t") {
			tcount = atoi(argv[++i]);
			if(tcount < 1)
				tcount = 1;
		} else if(arg == "-m") {
			method = argv[++i];
		} else {
			if(i < argc - 1) {
				files.push_back(arg);
				bands.push_back(atoi(argv[++i]));
			} else {
				outfile = arg;
			}
		}
	}

	std::string targetFile = files.back();
	int targetBand = bands.back();

	std::string filebase = outfile.substr(0, outfile.find('.'));
	std::cout << filebase << "\n";

	Ctx rdiff;

	{

		Ctx data1;
		loadRaster(targetFile, targetBand, data1);

		Ctx data2;
		loadRasters(files, bands, files.size() - 1, 1, data2);

		bool hasMask = !maskfile.empty();
		Ctx mask;
		if(hasMask)
			loadRaster(maskfile, maskband, mask);

		bool tonan;
		for(int row1 = 0; row1 < data1.rows; ++row1) {
			for(int col1 = 0; col1 < data1.cols; ++col1) {
				tonan = true;
				double v1 = data1.data[row1 * data1.cols + col1];
				if(v1 != data1.nd) {
					double x = toX(col1, data1.trans);
					double y = toY(row1, data1.trans);
					int col2 = toCol(x, data2.trans);
					int row2 = toRow(y, data2.trans);
					int mc = toCol(x, mask.trans);
					int mr = toRow(y, mask.trans);
					if(!(col2 < 0 || row2 < 0 || col2 >= data2.cols || row2 >= data2.rows) && !(mc < 0 || mr < 0 || mc >= mask.cols || mr >= mask.rows)) {
						if(!hasMask || mask.data[mr * mask.cols + mc] == 1) {
							double v2 = data2.data[row2 * data2.cols + col2];
							if(v2 != data2.nd) {
								data1.data[row1 * data1.cols + col1] = v2 - v1;
								tonan = false;
							}
						}
					}
				}
				if(tonan)
					data1.data[row1 * data1.cols + col1] = std::nan("");
			}
		}

		if(diff) {
			saveRaster(filebase + "_diff.tif", data1);
			saveRaster("/tmp/data2.tif", data2);
		}

		if(resample == 1) {
			rdiff = data1;
			rdiff.data.swap(data1.data);
		} else {
			rdiff = data1;
			rdiff.cols = (int) std::ceil((float) rdiff.cols / resample);
			rdiff.rows = (int) std::ceil((float) rdiff.rows / resample);
			rdiff.trans[1] *= resample;
			rdiff.trans[5] *= resample;
			rdiff.data.resize(rdiff.cols * rdiff.rows);
			std::fill(rdiff.data.begin(), rdiff.data.end(), 0);
			std::vector<int> counts(rdiff.cols * rdiff.rows);
			std::fill(counts.begin(), counts.end(), 0);
			float v;
			for(int row = 0; row < data1.rows; ++row) {
				for(int col = 0; col < data1.cols; ++col) {
					int cc = col / resample;
					int rr = row / resample;
					if(!std::isnan(v = data1.data[row * data1.cols + col])) {
						rdiff.data[rr * rdiff.cols + cc] += v;
						++counts[rr * rdiff.cols + cc];
					}
				}
			}
			for(size_t i = 0; i < rdiff.data.size(); ++i) {
				if(counts[i]) {
					rdiff.data[i] /= counts[i];
				} else {
					rdiff.data[i] = std::nan("");
				}
			}
		}
	}

	{

		Ctx dst = rdiff;
		std::fill(dst.data.begin(), dst.data.end(), 0);

		std::list<int> rowq;
		for(int row = 0; row < rdiff.rows; ++row)
			rowq.push_back(row);

		std::mutex qmtx;
		std::mutex dmtx;
		std::vector<std::thread> threads;

		int size = (radius * 2) / std::abs(rdiff.trans[1]);
		if(size % 2 == 0)
			size++;
		if(method == "gauss") {
			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processGauss, &rowq, &qmtx, &rdiff, &dst, &dmtx, size, (float) size * 0.25);
		} else if(method == "cosine") {
			float cos[1001];
			for(int i = 0; i <= 1000; ++i)
				cos[i] = std::cos((float) i / 1000 * M_PI) / 2.0 + 0.5;
			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processCos, &rowq, &qmtx, &rdiff, &dst, &dmtx, size, cos);
		} else if(method == "idw") {
			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processIDW, &rowq, &qmtx, &rdiff, &dst, &dmtx, size);
		} else if(method == "dw") {
			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processDW, &rowq, &qmtx, &rdiff, &dst, &dmtx, size);
		} else {
			std::cerr << "Unknown method: " << method << "\n";
			return 1;
		}

		for(int i = 0; i < tcount; ++i) {
			if(threads[i].joinable())
				threads[i].join();
		}

		rdiff.data.swap(dst.data);
	}

	Ctx final;
	loadRaster(targetFile, targetBand, final);

	double v;
	for(int row = 0; row < final.rows; ++row) {
		for(int col = 0; col < final.cols; ++col) {
			if((v = final.data[row * final.cols + col]) != final.nd) {
				float x = toX(col, final.trans);
				float y = toY(row, final.trans);
				int c = toCol(x, rdiff.trans);
				int r = toRow(y, rdiff.trans);
				float f = rdiff.data[r * rdiff.cols + c];
				final.data[row * final.cols + col] = v + f;
			}
		}
	}

	if(diff)
		saveRaster(filebase + "_diff_smooth.tif", rdiff);

	saveRaster(outfile, final);

	return 0;
}

