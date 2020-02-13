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
#include <unordered_set>

#include <gdal_priv.h>

#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>

class Ctx {
private:
	std::vector<float> _data;

public:
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
		_data.resize(other.size());
		for(int i = 0; i < 6; ++i)
			trans[i] = other.trans[i];
		return *this;
	}

	size_t size() const {
		return _data.size();
	}

	void resize(size_t i) {
		_data.resize(i);
	}

	std::vector<float>& data() {
		return _data;
	}

	float& get(size_t i) {
		if(i < _data.size())
			return _data[i];
		throw std::runtime_error("Invalid index.");
	}

	float minx() const {
		return trans[1] < 0 ? trans[0] + cols * trans[1] : trans[0];
	}

	float maxx() const {
		return trans[1] > 0 ? trans[0] + cols * trans[1] : trans[0];
	}

	float miny() const {
		return trans[5] < 0 ? trans[3] + rows * trans[5] : trans[3];
	}

	float maxy() const {
		return trans[5] > 0 ? trans[3] + cols * trans[5] : trans[3];
	}
};

int toCol(float x, double* trans) {
	return (int) (x - trans[0]) / trans[1];
}

int toRow(float y, double* trans) {
	return (int) (y - trans[3]) / trans[5];
}

float toX(int col, double* trans) {
	return (col * trans[1]) + trans[0] + trans[1] * 0.5;
}

float toY(int row, double* trans) {
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
	data.resize(data.cols * data.rows);
	data.projection = ds->GetProjectionRef();
	if(CE_None != bnd->RasterIO(GF_Read, 0, 0, data.cols, data.rows, data.data().data(), data.cols, data.rows, GDT_Float32, 0, 0, 0)) {
		GDALClose(ds);
		return false;
	}
	GDALClose(ds);
	return true;
}

bool loadRasters(const std::vector<std::string>& files, std::vector<int> bands, int count, int resample, Ctx& data) {

	GDALAllRegister();

	double trans[6];
	float minx = 99999999.0, miny = 999999999.0;
	float maxx = -99999999.0, maxy = -99999999.0;
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
	data.resize(data.cols * data.rows);
	data.nd = std::nan("");
	std::fill(data.data().begin(), data.data().end(), data.nd);
	for(int i = 0; i < 6; ++i)
		data.trans[i] = trans[i];

	std::vector<float> buf;
	std::vector<int> counts;

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
					float v = buf[r * (cols / resample) + c];
					if(v != nd)
						data.get(rr * data.cols + cc) = v;
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
	if(CE_None != ds->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, data.cols, data.rows, data.data().data(), data.cols, data.rows, GDT_Float32, 0, 0, 0)) {
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
	float v;
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

void processIDW(std::list<int>* rowq, std::mutex* qmtx, pcl::KdTreeFLANN<pcl::PointXYZ>* tree,
		Ctx* data, std::mutex* dmtx, int ncount) {
	int row;
	float ex = 2.0;
	while(!rowq->empty()) {
		{
			std::lock_guard<std::mutex> lk(*qmtx);
			if(!rowq->empty()) {
				row = rowq->front();
				rowq->pop_front();
				std::cout << "Row " << row << " of " << data->rows << "\n";
			} else {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
				continue;
			}
		}
		size_t count;
		std::vector<int> indices;
		std::vector<float> dist;
		for(int col = 0; col < data->cols; ++col) {
			float& v = data->get(row * data->cols + col);
			if(v != data->nd && !std::isnan(v)) {
				float x = toX(col, data->trans);
				float y = toY(row, data->trans);
				pcl::PointXYZ q(x, y, 0);
				if((count = tree->nearestKSearch(q, ncount, indices, dist))) {
					float s = 0;
					float w = 0;
					for(int idx : indices) {
						pcl::PointXYZ pt = tree->getInputCloud()->at(idx);
						float d0 = std::sqrt(std::pow(pt.x - x, 2.0) + std::pow(pt.y - y, 2.0));
						if(d0 == 0) {
							s = pt.z;
							w = 1;
							break;
						} else {
							float w0 = 1.0 / std::pow(d0, ex);
							s += pt.z * w0;
							w += w0;
						}
					}
					if(w) {
						std::lock_guard<std::mutex> lk(*dmtx);
						v += s / w;
					}
				}
			}
		}
	}
}

void processDW(std::list<int>* rowq, std::mutex* qmtx, pcl::KdTreeFLANN<pcl::PointXYZ>* tree,
		Ctx* data, std::mutex* dmtx, int ncount) {
	int row;
	while(!rowq->empty()) {
		{
			std::lock_guard<std::mutex> lk(*qmtx);
			if(!rowq->empty()) {
				row = rowq->front();
				rowq->pop_front();
				std::cout << "Row " << row << " of " << data->rows << "\n";
			} else {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
				continue;
			}
		}
		size_t count;
		std::vector<int> indices;
		std::vector<float> dist;
		for(int col = 0; col < data->cols; ++col) {
			float& v = data->get(row * data->cols + col);
			if(v != data->nd && !std::isnan(v)) {
				float x = toX(col, data->trans);
				float y = toY(row, data->trans);
				pcl::PointXYZ q(x, y, 0);
				if((count = tree->nearestKSearch(q, ncount, indices, dist))) {
					float maxd = 0;
					for(float d : dist)
						maxd = std::max(maxd, d);
					maxd = std::sqrt(maxd);
					float s = 0;
					float w = 0;
					for(int idx : indices) {
						pcl::PointXYZ pt = tree->getInputCloud()->at(idx);
						float d0 = std::sqrt(std::pow(pt.x - x, 2.0) + std::pow(pt.y - y, 2.0));
						float w0 = 1.0 - d0 / maxd;
						s += pt.z * w0;
						w += w0;
					}
					if(w) {
						std::lock_guard<std::mutex> lk(*dmtx);
						v += s / w;
					}
				}
			}
		}
	}
}

float median(std::vector<float>& vals) {
	if(!vals.empty()) {
		std::sort(vals.begin(), vals.end());
		size_t c = vals.size();
		if(c % 2 == 0) {
			return (vals[c / 2] + vals[c / 2 + 1]) / 2.0;
		} else {
			return vals[c / 2];
		}
	}
	return std::nan("");
}

size_t medstats(Ctx& data, int col, int row, int c, int r, float& dist, float& med) {
	int rad = 3;
	float v;
	std::vector<float> vals;
	for(int rr = r - rad; rr < r + rad + 1; ++rr) {
		for(int cc = c - rad; cc < c + rad + 1; ++cc) {
			if(cc < 0 || rr < 0 || cc >= data.cols || rr >= data.rows)
				continue;
			v = data.get(rr * data.cols + cc);
			if(v != data.nd && !std::isnan(v))
				vals.push_back(v);
		}
	}
	med = median(vals);
	dist = std::pow(c - col, 2) + std::pow(r - row, 2);
	return vals.size();
}

/*
void processMedian(std::list<int>* rowq, std::mutex* qmtx, Ctx* data, Ctx* med, std::mutex* dmtx, int size) {
	std::vector<float> rowBuf(data->cols);
	int row;
	while(!rowq->empty()) {
		{
			std::lock_guard<std::mutex> lk(*qmtx);
			if(!rowq->empty()) {
				row = rowq->front();
				rowq->pop_front();
				std::cout << "Row " << row << " of " << data->rows << "\n";
			} else {
				std::this_thread::sleep_for(std::chrono::milliseconds(1));
				continue;
			}
		}

		std::vector<float> values;
		values.reserve(std::pow(size * 2 + 1, 2));
		size_t idx = 0;
		float v;
		for(int col = 0; col < data->cols; ++col) {
			int(int r = row - size; r < row + size + 1; ++r) {
				for(int c = col - size; c < col + size + 1; ++c) {
					if(c < 0 || r < 0 || c >= data->cols || r >= data->rows ||
							mask->data[r * mask->rows + c] != 1) {
						v = data->data[r * data->cols + c];
						if(v != data->nd || !std::isnan(v)) {
							values[idx++] = v;
						}
					}
				}
			}
		}
		std::lock_guard<std::mutex> lk(*dmtx);
		std::memcpy(data->data.data() + row * data->cols, rowBuf.data(), data->cols);
	}
}
*/

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
					if(!std::isnan((v0 = src->get(rr * src->cols + cc)))) {
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
				dst->get(row * src->cols + col) = s / w;
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
					if(!std::isnan((v0 = src->get(rr * src->cols + cc)))) {
						float w0 = std::exp(-0.5 * (c * c + r * r) / (sigma * sigma));
						//std::cout << d << " " << c << " " << r << " " << (c * c + r * r) << " " << rad2 << "\n";
						s += v0 * w0;
						w += gain;
					}
				}
			}
			if(w > 0) {
				std::lock_guard<std::mutex> lk(*dmtx);
				dst->get(row * src->cols + col) = s / w;
			}
		}
	}
}

int main(int argc, char** argv) {

	if(argc < 6) {
		std::cerr << "Usage: rastermerge [options] <<anchor file 1> <anchor band 1> [<anchor file 2> <anchor band 2> [...]]> <target file 2> <target band 2> <output file>\n"
				<< " -s <size>          The size of the median window in pixels.\n"
				<< " -c <count>         Number of neighbours to consider.\n"
				<< " -t <threads>       The number of threads.\n"
				<< " -m <method>        The method: idw, dw, gauss, cosine. Default IDW.\n"
				<< " -k <mask> <band>   A mask file. Pixel value 1 is kept.\n";
		return 1;
	}

	std::vector<std::string> files;
	std::vector<int> bands;
	int size = 11;
	int count = 20;
	int tcount = 4;
	std::string method = "dw";
	std::string outfile;
	std::string maskfile;
	int maskband = 0;

	for(int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-s") {
			size = atoi(argv[++i]);
			if(size % 2 == 0)
				++size;
		} else if(arg == "-k") {
			maskfile = argv[++i];
			maskband = atoi(argv[++i]);
		} else if(arg == "-c") {
			count = atoi(argv[++i]);
			if(count < 1)
				count = 20;
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

	Ctx data1;

	pcl::KdTreeFLANN<pcl::PointXYZ> tree;
	pcl::PointCloud<pcl::PointXYZ>::Ptr pts(new pcl::PointCloud<pcl::PointXYZ>);

	{

		loadRaster(targetFile, targetBand, data1);

		Ctx data2;
		loadRasters(files, bands, files.size() - 1, 1, data2);

		bool hasMask = !maskfile.empty();
		Ctx mask;
		if(hasMask)
			loadRaster(maskfile, maskband, mask);

		int cols = (int) std::ceil((float) data1.cols / size);
		int rows = (int) std::ceil((float) data1.rows / size);
		pts->width = cols;
		pts->height = rows;
		pts->points.resize(pts->width * pts->height);

		std::vector<float> values1;
		std::vector<float> values2;

		for(int row1 = 0; row1 < data1.rows; row1 += size) {
			for(int col1 = 0; col1 < data1.cols; col1 += size) {

				for(int r1 = row1 - size / 2; r1 < row1 + size / 2 + 1; ++r1) {
					for(int c1 = col1 - size / 2; c1 < col1 + size / 2 + 1; ++c1) {

						if(!(c1 < 0 || r1 < 0 || c1 >= data1.cols || r1 >= data1.rows)) {

							float v1 = data1.get(r1 * data1.cols + c1);
							if(v1 != data1.nd && !std::isnan(v1)) {

								float x = toX(c1, data1.trans);
								float y = toY(r1, data1.trans);
								int c2 = toCol(x, data2.trans);
								int r2 = toRow(y, data2.trans);
								int mc = toCol(x, mask.trans);
								int mr = toRow(y, mask.trans);

								if(!(c2 < 0 || r2 < 0 || c2 >= data2.cols || r2 >= data2.rows)
										&& (!hasMask || !(mc < 0 || mr < 0 || mc >= mask.cols || mr >= mask.rows))) {

									if(!hasMask || mask.get(mr * mask.cols + mc) == 1) {

										float v2 = data2.get(r2 * data2.cols + c2);
										if(v2 != data2.nd && !std::isnan(v2)) {
											values1.push_back(v1);
											values2.push_back(v2);
										}
									}
								}
							}
						}
					}
				}
				if(!values1.empty()) {
					float x = toX(col1, data1.trans);
					float y = toY(row1, data1.trans);
					float m1 = median(values1);
					float m2 = median(values2);
					pts->points.emplace_back(x, y, m2 - m1);
					values1.clear();
					values2.clear();
				}
			}
		}

		tree.setInputCloud(pts);
	}

	{

		std::list<int> rowq;
		for(int row = 0; row < data1.rows; ++row)
			rowq.push_back(row);

		std::mutex qmtx;
		std::mutex dmtx;
		std::vector<std::thread> threads;

		if(method == "gauss") {
//			for(int i = 0; i < tcount; ++i)
//				threads.emplace_back(processGauss, &rowq, &qmtx, &tree, &dst, &dmtx, size, (float) size * 0.25);
		} else if(method == "cosine") {
			/*
			float cos[1001];
			for(int i = 0; i <= 1000; ++i)
				cos[i] = std::cos((float) i / 1000 * M_PI) / 2.0 + 0.5;
			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processCos, &rowq, &qmtx, &tree, &data1, &dmtx, size, cos);
				*/
		} else if(method == "idw") {

			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processIDW, &rowq, &qmtx, &tree, &data1, &dmtx, count);

		} else if(method == "dw") {

			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processDW, &rowq, &qmtx, &tree, &data1, &dmtx, count);

		} else if(method == "med") {

			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processDW, &rowq, &qmtx, &tree, &data1, &dmtx, count);

		} else {
			std::cerr << "Unknown method: " << method << "\n";
			return 1;
		}

		for(std::thread& t : threads) {
			if(t.joinable())
				t.join();
		}
	}

	saveRaster(outfile, data1);

	return 0;
}

