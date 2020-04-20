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
#include <fstream>

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include "util.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3   Point;
typedef Delaunay::Face_handle Face_handle;
typedef Delaunay::Vertex Vertex;

/**
 * Contains a raster grid and related data/methods.
 */
class Ctx {
private:
	std::vector<float> _data;

public:
	int cols;
	int rows;
	double trans[6];
	float nd;
	std::string projection;
/*
	static Ctx fromTriangulation(const std::string& ptsfile, const std::string& layer,
			const std::string& idcol, const std::vector<int>& ids,
			float resx, float resy, const std::string& prprojection) {

		GDALAllRegister();

		GDALDataset* ds = (GDALDataset*) GDALOpenEx(ptsfile.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL);
		if(!ds)
			g_runerr("Couldn't open database.");

		std::stringstream ss;
		ss << "SELECT * FROM " << layer << " WHERE " << idcol << " IN (";
		for(size_t i = 0; i < ids.size(); ++i) {
			if(i > 0)
				ss << ","
			ss << ids[i];
		}
		ss << ")";

		OGRLayer* lyr = ds->ExecuteSQL(ss.str().c_str(), NULL, "SQLite");
		if(!lyr)
			g_runerr("Couldn't get layer: " << layer);

		std::vector<Point> points;
		for(OGRFeature& feat : lyr) {
			OGRGeometry* g = feat.GetGeometryRef();
			OGRPoint* p = g->toPoint();
			data.get(p->getX(), p->getY());
		}

		ds->ReleaseResultSet(lyr);
		GDALClose(ds);
		  for(int r = 0; r < rows; ++r) {
			  for(int c = 0; c < cols; ++c) {
				  float& v = get(c, r);
				  if(v != nd && !std::isnan(v))
					  points.emplace_back(toX(c), toY(r), v);
			  }
		  }

		  int rad = std::max(cols, rows) * 100;
		  for(int i = 0; i < 360; ++i) {
			  int c = cols / 2 + std::cos((float) i * M_PI / 180.) * rad;
			  int r = rows / 2 + std::sin((float) i * M_PI / 180.) * rad;
			  points.emplace_back(toX(c), toY(r), 0);
		  }

		  Delaunay d(points.begin(), points.end());
		  Face_handle h;
		  std::vector<float> vertices(9);

		  for(int r = 0; r < rows; ++r) {
			  for(int c = 0; c < cols; ++c) {
				  float& v = get(c, r);
				  if(v == nd || std::isnan(v)) {
					  float x = toX(c);
					  float y = toY(r);
					  h = d.locate(Point(x, y, 0), h);
					  for(int i = 0; i < 3; ++i) {
						  Point t = h->vertex(i)->point();
						  vertices[i * 3 + 0] = t.x();
						  vertices[i * 3 + 1] = t.y();
						  vertices[i * 3 + 2] = t.z();
					  }
					  v = bary(x, y, vertices);
				  }
			  }
		  }
	}
*/

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

	Ctx& operator+=(float v) {
		for(float& f : data()) {
			if(f != nd && !std::isnan(f))
				f += v;
		}
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

	void voidfill() {
		int ct;
		do {
			ct = 0;
			for(int row = 0; row < rows; ++row) {
				for(int col = 0; col < cols; ++col) {
					float& v = get(col, row);
					if(v == nd) {
						int c = 0;
						float s = 0;
						for(int r = row - 1; r < row + 2; ++r) {
							for(int c = col - 1; c < col + 2; ++c) {
								if(c >= 0 && r >= 0 && c < cols && r < rows) {
									float& vv = get(c, r);
									if(vv != nd) {
										s += vv;
										++c;
									}
								}
							}
						}
						if(c) {
							v = s / c;
							++ct;
						}
					}
				}
			}
		} while(ct);
	}

	float bary(float x, float y, const std::vector<float>& v) {
		float x1 = v[0];
		float y1 = v[1];
		float z1 = v[2];
		float x2 = v[3];
		float y2 = v[4];
		float z2 = v[5];
		float x3 = v[6];
		float y3 = v[7];
		float z3 = v[8];
		float w1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
		float w2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
		float w3 = 1 - w1 - w2;
		float z = (z1 * w1) + (z2 * w2) + (z3 * w3);
		return z;
	}

	void voidfilldelaunay() {
	  std::vector<Point> points;
	  for(int r = 0; r < rows; ++r) {
		  for(int c = 0; c < cols; ++c) {
			  float& v = get(c, r);
			  if(v != nd && !std::isnan(v))
				  points.emplace_back(toX(c), toY(r), v);
		  }
	  }

	  int rad = std::max(cols, rows) * 100;
	  for(int i = 0; i < 360; ++i) {
		  int c = cols / 2 + std::cos((float) i * M_PI / 180.) * rad;
		  int r = rows / 2 + std::sin((float) i * M_PI / 180.) * rad;
		  points.emplace_back(toX(c), toY(r), 0);
	  }

	  Delaunay d(points.begin(), points.end());
	  Face_handle h;
	  std::vector<float> vertices(9);

	  for(int r = 0; r < rows; ++r) {
		  for(int c = 0; c < cols; ++c) {
			  float& v = get(c, r);
			  if(v == nd || std::isnan(v)) {
				  float x = toX(c);
				  float y = toY(r);
				  h = d.locate(Point(x, y, 0), h);
				  for(int i = 0; i < 3; ++i) {
					  Point t = h->vertex(i)->point();
					  vertices[i * 3 + 0] = t.x();
					  vertices[i * 3 + 1] = t.y();
					  vertices[i * 3 + 2] = t.z();
				  }
				  v = bary(x, y, vertices);
			  }
		  }
	  }
	}

	void resample(float scale) {
		cols = (int) std::ceil(cols / scale);
		rows = (int) std::ceil(rows / scale);
		trans[1] *= scale;
		trans[5] *= scale;
		_data.resize(cols * rows);
	}

	float& get(size_t i) {
		if(i < _data.size())
			return _data[i];
		throw std::runtime_error("Invalid index.");
	}

	float& get(int col, int row) {
		return get(row * cols + col);
	}

	float& get(float x, float y) {
		return get(toCol(x), toRow(y));
	}

	float minx() const {
		return trans[1] < 0 ? trans[0] + (cols + 1) * trans[1] : trans[0];
	}

	float maxx() const {
		return trans[1] > 0 ? trans[0] + (cols + 1) * trans[1] : trans[0];
	}

	float miny() const {
		return trans[5] < 0 ? trans[3] + (rows + 1) * trans[5] : trans[3];
	}

	float maxy() const {
		return trans[5] > 0 ? trans[3] + (rows + 1) * trans[5] : trans[3];
	}

	int toCol(float x) {
		return (int) (x - trans[0]) / trans[1];
	}

	int toRow(float y) {
		return (int) (y - trans[3]) / trans[5];
	}

	float toX(int col) {
		return (col * trans[1]) + trans[0] + trans[1] * 0.5;
	}

	float toY(int row) {
		return (row * trans[5]) + trans[3] + trans[5] * 0.5;
	}

	int nextCol(int col) {
		if(trans[1] > 0) {
			if(++col >= cols)
				return -1;
		} else if(--col < 0) {
			return -1;
		}
		return col;
	}

	int prevCol(int col) {
		if(trans[1] > 0) {
			if(--col < 0)
				return -1;
		} else if(++col >= cols) {
			return -1;
		}
		return col;
	}

	int nextRow(int row) {
		if(trans[5] > 0) {
			if(++row >= rows)
				return -1;
		} else if(--row < 0) {
			return -1;
		}
		return row;
	}

	int prevRow(int row) {
		if(trans[5] > 0) {
			if(--row < 0)
				return -1;
		} else if(++row >= rows) {
			return -1;
		}
		return row;
	}

	bool load(const std::string& file, int band) {

		GDALAllRegister();

		GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(file.c_str(), GA_ReadOnly));
		if(!ds)
			return false;
		ds->GetGeoTransform(trans);
		cols = ds->GetRasterXSize();
		rows = ds->GetRasterYSize();
		ds->GetGeoTransform(trans);
		GDALRasterBand* bnd = ds->GetRasterBand(band);
		nd = bnd->GetNoDataValue();
		resize(cols * rows);
		projection = ds->GetProjectionRef();
		if(CE_None != bnd->RasterIO(GF_Read, 0, 0, cols, rows, data().data(), cols, rows, GDT_Float32, 0, 0, 0)) {
			GDALClose(ds);
			return false;
		}
		GDALClose(ds);
		return true;
	}

	bool load(const std::vector<std::string>& files, std::vector<int> bands, int count, int resample) {

		GDALAllRegister();

		float minx, miny, maxx, maxy;

		minx = miny = G_DBL_MAX_POS;
		maxx = maxy = G_DBL_MAX_NEG;

		for(size_t i = 0; i < (size_t) count; ++i) {
			std::string file = files[i];
			int band = bands[i];

			GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(file.c_str(), GA_ReadOnly));
			if(!ds)
				return false;

			double trans0[6];
			ds->GetGeoTransform(trans0);

			int cs = ds->GetRasterXSize();
			int rs = ds->GetRasterYSize();

			if(i == 0) {
				nd = ds->GetRasterBand(band)->GetNoDataValue();
				projection = ds->GetProjectionRef();
				for(int j = 0; j < 6; ++j)
					trans[j] = trans0[j];
			}

			if(trans0[1] > 0) {
				trans[0] = std::min(trans[0], trans0[0]);
				minx = std::min(minx, (float) trans0[0]);
				maxx = std::max(maxx, (float) (trans0[0] + cs * trans0[1]));
			} else {
				trans[0] = std::max(trans[0], trans0[0]);
				maxx = std::max(maxx, (float) trans0[0]);
				minx = std::min(minx, (float) (trans0[0] + cs * trans0[1]));
			}
			if(trans0[5] > 0) {
				trans[3] = std::min(trans[3], trans0[3]);
				miny = std::min(miny, (float) trans0[3]);
				maxy = std::max(maxy, (float) (trans0[3] + rs * trans0[5]));
			} else {
				trans[3] = std::max(trans[3], trans0[3]);
				maxy = std::max(maxy, (float) trans0[3]);
				miny = std::min(miny, (float) (trans0[3] + rs * trans0[5]));
			}
			GDALClose(ds);
		}

		trans[1] *= resample;
		trans[5] *= resample;
		cols = (int) std::ceil((maxx - minx) / std::abs(trans[1]));
		rows = (int) std::ceil((maxy - miny) / std::abs(trans[5]));

		resize(cols * rows);
		std::fill(_data.begin(), _data.end(), 0);

		std::vector<float> buf;
		std::vector<int> counts(cols * rows);
		std::fill(counts.begin(), counts.end(), 0);

		for(size_t i = 0; i < (size_t) count; ++i) {
			std::string file = files[i];
			int band = bands[i];

			GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(file.c_str(), GA_ReadOnly));
			if(!ds)
				return false;

			int cs = ds->GetRasterXSize();
			int rs = ds->GetRasterYSize();

			double trans0[6];
			ds->GetGeoTransform(trans0);
			trans0[1] *= resample;
			trans0[5] *= resample;

			int csr = (int) std::ceil((float) cs / resample);
			int rsr = (int) std::ceil((float) rs / resample);

			buf.resize(csr * rsr);

			GDALRasterBand* bnd = ds->GetRasterBand(band);
			if(CE_None != bnd->RasterIO(GF_Read, 0, 0, cs, rs, buf.data(), csr, rsr, GDT_Float32, 0, 0, 0)) {
				GDALClose(ds);
				return false;
			}

			float nd0 = bnd->GetNoDataValue();
			int co = toCol(trans0[0]);
			int ro = toRow(trans0[3]);
			for(int r = 0; r < rsr; ++r) {
				for(int c = 0; c < csr; ++c) {
					float x = toX(c + co);
					float y = toY(r + ro);
					int cc = toCol(x);
					int rr = toRow(y);
					if(!(cc < 0 || rr < 0 || cc >= cols || rr >= rows)) {
						float v = buf[r * csr + c];
						if(v != nd0) {
							float& vv = get(cc, rr);
							vv += v;
							counts[rr * cols + cc]++;
						}
					}
				}
			}

			GDALClose(ds);
		}

		int ct;
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				float& v = get(c, r);
				if((ct = counts[r * cols + c]) > 0) {
					v /= ct;
				} else {
					v = nd;
				}
			}
		}

		return true;
	}

	bool save(const std::string& file) {
		GDALDriverManager* dm = GetGDALDriverManager();
		GDALDriver* drv = dm->GetDriverByName("GTiff");
		GDALDataset* ds = drv->Create(file.c_str(), cols, rows, 1, GDT_Float32, 0);
		if(!ds)
			return false;
		ds->SetProjection(projection.c_str());
		ds->SetGeoTransform(trans);
		ds->GetRasterBand(1)->SetNoDataValue(nd);
		if(CE_None != ds->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, cols, rows, data().data(), cols, rows, GDT_Float32, 0, 0, 0)) {
			GDALClose(ds);
			return false;
		}
		GDALClose(ds);
		return true;
	}

	bool isEdgePixel(int col, int row, int size, float prop = 0.30) {
		float v;
		if(col < 0 || row < 0 || col >= cols || row >= rows || (v = get(col, row)) == nd || std::isnan(v))
			return false;
		int ct = 0, t = 0;
		for(int r = row - size; r < row + size + 1; ++r) {
			for(int c = col - 1; c < col + 2; ++c) {
				if(c < 0 || r < 0 || c >= cols || r >= rows)
					return true;
				++t;
				if((v = get(c, r)) == nd || std::isnan(v))
					++ct;
			}
		}
		return (float) ct / t >= prop;
	}

	// Get the median for the cell neighbourhood with radius.
	// Return nan if something goes wrong.
	float median(int col, int row, float radius) {
		if(col < 0 || row < 0 || col >= cols || row >= rows)
			return std::nan("");
		int o = std::max(1, (int) std::ceil(radius / std::abs(trans[1])));
		static std::vector<float> values;
		values.resize(0);
		float v;
		for(int r = row - o; r < row + o + 1; ++r) {
			for(int c = col - o; c < col + o + 1; ++c) {
				if(!(c < 0 || r < 0 || c >= cols || r >= rows)
						&& (v = get(c, r)) != nd && !std::isnan(v)) {
					values.push_back(v);
				}
			}
		}
		if(values.size() < 2) {
			return std::nan("");
		} else {
			std::sort(values.begin(), values.end());
			if(values.size() % 2 == 0) {
				return (values[values.size() / 2 - 1] + values[values.size() / 2]) / 2.0;
			} else {
				return values[values.size() / 2];
			}
		}
	}

	void fill(float f) {
		std::fill(_data.begin(), _data.end(), f);
	}
};





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

void processDW(std::list<int>* rowq, bool idw, float exponent, std::mutex* qmtx, pcl::KdTreeFLANN<pcl::PointXYZ>* tree,
		Ctx* data, Ctx* diff, std::mutex* dmtx, float radius, int ncount) {

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
				float x = data->toX(col);
				float y = data->toY(row);
				pcl::PointXYZ q(x, y, 0);
				if((count = tree->radiusSearch(q, radius, indices, dist, 0)) > (size_t) ncount) {
					float s = 0;
					float w = 0;
					for(size_t i = 0; i < indices.size(); ++i) {
						size_t idx = indices[i];
						pcl::PointXYZ pt = tree->getInputCloud()->at(idx);
						float d0 = dist[i];
						if(idw) {
							if(d0 == 0) {
								s = pt.z;
								w = 1;
								break;
							} else {
								float w0 = 1.0 / (exponent == 2 ? d0 : std::pow(std::sqrt(d0), exponent));
								s += pt.z * w0;
								w += w0;
							}
						} else {
							float d = std::sqrt(dist[i]);
							float w0 = 1.0 - std::min(d / radius, 1.0f);
							s += pt.z * w0;
							w += 1; // TODO: CHooseablew0;
						}
					}
					if(w) {
						std::lock_guard<std::mutex> lk(*dmtx);
						v += s / w;
						diff->get(row * data->cols + col) = s / w;
					}
				}
			}
		}
	}
}

void processBI(std::list<int>* rowq, std::mutex* qmtx, Ctx* medgrid,
		Ctx* data, Ctx* diff, std::mutex* dmtx) {
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
		std::vector<int> indices;
		std::vector<float> dist;
		for(int col = 0; col < data->cols; ++col) {
			float& v = data->get(row * data->cols + col);
			if(v != data->nd && !std::isnan(v)) {
				float x = data->toX(col);
				float y = data->toY(row);
				int mc = medgrid->toCol(x);
				int mr = medgrid->toRow(y);
				float x1 = medgrid->toX(mc);
				float y1 = medgrid->toY(mr);
				float x2, y2;
				int mc0, mr0;
				if(x > x1) {
					mc0 = medgrid->nextCol(mc);
				} else {
					mc0 = medgrid->prevCol(mc);
				}
				if(y > y1) {
					mr0 = medgrid->nextRow(mr);
				} else {
					mr0 = medgrid->prevRow(mr);
				}
				if(mr0 < 0 || mc0 < 0)
					continue;

				x1 = medgrid->toX(mc), x2 = medgrid->toX(mc0);
				y1 = medgrid->toY(mr), y2 = medgrid->toY(mr0);

				float z1 = medgrid->get(x1, y1);
				float z2 = medgrid->get(x2, y1);
				float z3 = medgrid->get(x1, y2);
				float z4 = medgrid->get(x2, y2);

				float z11 = z1 + (z2 - z1) * ((x - x1) / (x2 - x1));
				float z22 = z3 + (z4 - z3) * ((x - x1) / (x2 - x1));
				float z = z11 + (z22 - z11) * ((y - y1) / (y2 - y1));
				std::lock_guard<std::mutex> lk(*dmtx);
				if(x == 475835 && y == 6500105)
					std::cout << z << "\n";
				diff->get(x, y) = z;
				v += z;
			}
		}
	}
}


void processSpline(geo::util::BivariateSpline* spline, std::list<int>* rowq,
		std::mutex* qmtx, Ctx* data, std::mutex* dmtx) {
	/*
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
		for(int col = 0; col < data->cols; ++col) {
			float& v = data->get(row * data->cols + col);
			if(v != data->nd && !std::isnan(v)) {
				float x = data->toX(col);
				float y = data->toY(row);
				float z[2];
				if(!spline->evaluate(&x, 1, &y, 1, z, 2)) {
					std::lock_guard<std::mutex> lk(*dmtx);
					v += z[0];
				}
			}
		}
	}
	*/
}

float median(std::vector<float>& vals) {
	if(vals.size() >= 2) {
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

float mean(std::vector<float>& vals) {
	if(!vals.empty()) {
		float x = 0;
		for(float v : vals)
			x += v;
		return x / vals.size();
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


void processGauss(std::list<int>* rowq, std::mutex* qmtx, pcl::KdTreeFLANN<pcl::PointXYZ>* tree, Ctx* data, Ctx* diff, std::mutex* dmtx, float sigma) {
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

		float rad = sigma * std::abs(data->trans[1]) * 3;
		int count;
		std::vector<float> dist;
		std::vector<int> idx;
		for(int col = 0; col < data->cols; ++col) {
			float& v = data->get(col, row);
			if(v == data->nd || std::isnan(v))
				continue;
			float x = data->toX(col);
			float y = data->toY(row);
			pcl::PointXYZ q(x, y , 0);
			if((count = tree->radiusSearch(q, rad, idx, dist, (unsigned int) 0))) {
				float df = 0;
				for(size_t i = 0; i < idx.size(); ++i) {
					const pcl::PointXYZ& pt = tree->getInputCloud()->at(idx[i]);
					float w0 = std::exp(-0.5 * (float) dist[i] / (sigma * sigma));
					df += pt.z * w0;
				}
				std::lock_guard<std::mutex> lk(*dmtx);
				diff->get(col, row) = df;
			}
		}
	}
}

int main(int argc, char** argv) {

	if(argc < 6) {
		std::cerr << "Usage: rastermerge [options] <<anchor file 1> <anchor band 1> [<anchor file 2> <anchor band 2> [...]]> <target file 2> <target band 2> <output file>\n"
				<< " -s <size>          The step size in pixels.\n"
				<< " -c <count>         Number of neighbours to consider.\n"
				<< " -r <radius>        Search radius.\n"
				<< " -d <radius>        Radius for median calculation.\n"
				<< " -t <threads>       The number of threads.\n"
				<< " -m <method>        The method: idw, dw, gauss, cosine. Default IDW.\n"
				<< " -e <exponent>      Exponent for IDW (default 2).\n"
				<< " -k <mask> <band>   A mask file. Pixel value 1 is kept.\n";
		return 1;
	}

	std::vector<std::string> files;
	std::vector<int> bands;
	int size = 11;
	int count = 0;
	int tcount = 4;
	std::string method = "dw";
	std::string outfile;
	std::string maskfile;
	int maskband = 0;
	float radius = 100;
	float mradius = 100;
	float exponent = 2;
	bool edges = false;

	for(int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-s") {
			size = atoi(argv[++i]);
		} else if(arg == "-k") {
			maskfile = argv[++i];
			maskband = atoi(argv[++i]);
		} else if(arg == "-c") {
			count = atoi(argv[++i]);
		} else if(arg == "-x") {
			edges = true;
		} else if(arg == "-e") {
			exponent = atof(argv[++i]);
		} else if(arg == "-t") {
			tcount = atoi(argv[++i]);
			if(tcount < 1)
				tcount = 1;
		} else if(arg == "-m") {
			method = argv[++i];
		} else if(arg == "-r") {
			radius = atof(argv[++i]);
		} else if(arg == "-d") {
			mradius = atof(argv[++i]);
		} else {
			if(i < argc - 1) {
				files.push_back(arg);
				bands.push_back(atoi(argv[++i]));
			} else {
				outfile = arg;
			}
		}
	}

	//std::cout << "Radius: " << radius << ", size " << size << ", window: " << mradius << ", exponent " << exponent << ", count " << count << "\n";
	//return 1;

	std::string targetFile = files.back();
	int targetBand = bands.back();

	Ctx data1;
	Ctx medgrid;

	pcl::KdTreeFLANN<pcl::PointXYZ> tree;
	pcl::PointCloud<pcl::PointXYZ>::Ptr pts(new pcl::PointCloud<pcl::PointXYZ>);
//	geo::util::BivariateSpline bspline;
	std::ofstream ptsf("pts.csv");

	{

		data1.load(targetFile, targetBand);

		Ctx data2;
		data2.load(files, bands, files.size() - 1, 1);
		data2.save("/tmp/tmp_data2.tif");

		bool hasMask = !maskfile.empty();
		Ctx mask;
		if(hasMask)
			mask.load(maskfile, maskband);

		int cols = (int) std::ceil((float) data1.cols / size);
		int rows = (int) std::ceil((float) data1.rows / size);

		std::vector<float> values1;
		std::vector<float> values2;

		std::vector<float> bvx;
		std::vector<float> bvy;
		std::vector<float> bvz;
		std::vector<float> bvw;

		std::list<std::tuple<float, float, float>> medpts;

		medgrid = data1;
		medgrid.resample(size);
		medgrid.fill(medgrid.nd);

		if(method == "bi") {
		} else if(method == "dw" || method == "idw" || method == "gauss") {
			pts->width = cols;
			pts->height = rows;
			pts->points.resize(pts->width * pts->height);
		}

		float meandif = 0;

		{
			int meanct = 0;

			for(int r1 = 0; r1 < data1.rows; r1 += size) {
				for(int c1 = 0; c1 < data1.cols; c1 += size) {

					float v1 = data1.median(c1, r1, mradius);
					if(v1 != data1.nd && !std::isnan(v1)) {

						float x = data1.toX(c1);
						float y = data1.toY(r1);
						int c2 = data2.toCol(x);
						int r2 = data2.toRow(y);
						int mc = hasMask ? mask.toCol(x) : 0;
						int mr = hasMask ? mask.toRow(y) : 0;

						if(!(c2 < 0 || r2 < 0 || c2 >= data2.cols || r2 >= data2.rows)
								&& (!hasMask || !(mc < 0 || mr < 0 || mc >= mask.cols || mr >= mask.rows))) {

							if(hasMask && mask.get(mr * mask.cols + mc) != 1)
								continue;

							if(edges && !(data1.isEdgePixel(c1, r1, 5, 0) || data2.isEdgePixel(c2, r2, 5, 0)))
								continue;

							float v2 = data2.median(c2, r2, mradius);
							if(v2 != data2.nd && !std::isnan(v2)) {
								meandif += v2 - v1;
								++meanct;
							}
						}
					}
				}
			}
			data1 += meandif / meanct;
		}

		for(int r1 = 0; r1 < data1.rows; r1 += size) {
			for(int c1 = 0; c1 < data1.cols; c1 += size) {

				float v1 = data1.median(c1, r1, mradius);
				if(v1 != data1.nd && !std::isnan(v1)) {

					float x = data1.toX(c1);
					float y = data1.toY(r1);
					int c2 = data2.toCol(x);
					int r2 = data2.toRow(y);
					int mc = hasMask ? mask.toCol(x) : 0;
					int mr = hasMask ? mask.toRow(y) : 0;

					if(!(c2 < 0 || r2 < 0 || c2 >= data2.cols || r2 >= data2.rows)
							&& (!hasMask || !(mc < 0 || mr < 0 || mc >= mask.cols || mr >= mask.rows))) {

						if(hasMask && mask.get(mr * mask.cols + mc) != 1)
							continue;

						if(edges && !(data1.isEdgePixel(c1, r1, 5, 0) || data2.isEdgePixel(c2, r2, 5, 0)))
							continue;

						float v2 = data2.median(c2, r2, mradius);
						if(v2 != data2.nd && !std::isnan(v2)) {
							values1.push_back(v1);
							values2.push_back(v2);
						}

					}
				}

				float x = data1.toX(c1);
				float y = data1.toY(r1);
				std::string meanmed = "mean";
				if(!values1.empty()) {
					float m1 = meanmed == "mean" ? mean(values1) : median(values1);
					float m2 = meanmed == "mean" ? mean(values2) : median(values2);
					if(method == "spline") {
/*						bvx.push_back(x);
						bvy.push_back(y);
						bvz.push_back(m2 - m1);
						bvw.push_back(1);*/
					} else if(method == "idw" || method == "dw" || method == "gauss") {
						pts->points.emplace_back(x, y, m2 - m1);
					} else if(method == "bi") {
						medpts.emplace_back(x, y, m2 - m1);
					}
					medgrid.get(x, y) = m2 - m1;
					ptsf << std::setprecision(6) << std::fixed << x << "," << y << "," << (m2-m1) << "\n";
					values1.clear();
					values2.clear();
				} else if(method == "spline") {
					/*
					bvx.push_back(x);
					bvy.push_back(y);
					bvz.push_back(0);
					bvw.push_back(1);
					*/
				}
			}
		}

		if(method == "spline") {
			/*
			float smooth = 100;
			bspline.init(smooth, bvx, bvy, bvz, bvw,
					(float) data1.minx(), (float) data1.miny(),
					(float) data1.maxx(), (float) data1.maxy());
					*/
		} else if(method == "idw" || method == "dw" || method == "gauss") {
			tree.setInputCloud(pts);
		} else if(method == "bi") {
			float s = 0;
			int c = 0;
			for(const auto& m : medpts) {
				s += std::get<2>(m);
				++c;
			}
			float mean = s / c;
			s = 0;
			for(const auto& m : medpts) {
				s += std::pow(std::get<2>(m) - mean, 2.0);
			}
			s = std::sqrt(s) / 2;
			for(const auto& m : medpts) {
				float z = std::get<2>(m);
				if(z <= s && z >= -s)
					medgrid.get(std::get<0>(m), std::get<1>(m)) = z;
			}
		}
	}

	Ctx diff = data1;
	diff.fill(0);

	{

		std::list<int> rowq;
		for(int row = 0; row < data1.rows; ++row)
			rowq.push_back(row);

		std::mutex qmtx;
		std::mutex dmtx;
		std::vector<std::thread> threads;

		if(method == "gauss") {

			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processGauss, &rowq, &qmtx, &tree, &data1, &diff, &dmtx, radius);

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
				threads.emplace_back(processDW, &rowq, true, exponent, &qmtx, &tree, &data1, &diff, &dmtx, radius, count);

		} else if(method == "dw") {

			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processDW, &rowq, false, exponent, &qmtx, &tree, &data1, &diff, &dmtx, radius, count);

		} else if(method == "med") {

			//for(int i = 0; i < tcount; ++i)
				//threads.emplace_back(processMedian, &rowq, &qmtx, &tree, &data1, &dmtx, count);

		} else if(method == "spline") {


//			for(int i = 0; i < tcount; ++i)
//				threads.emplace_back(processSpline, &bspline, &rowq, &qmtx, &data1, &dmtx);

		} else if(method == "bi") {

			medgrid.save(outfile + ".bi1.tif");
			medgrid.voidfilldelaunay();

			for(int i = 0; i < tcount; ++i)
				threads.emplace_back(processBI, &rowq, &qmtx, &medgrid, &data1, &diff, &dmtx);

		} else {
			std::cerr << "Unknown method: " << method << "\n";
			return 1;
		}

		for(std::thread& t : threads) {
			if(t.joinable())
				t.join();
		}
	}

	data1.save(outfile);
	diff.save(outfile + ".diff.tif");
	medgrid.save(outfile + ".med.tif");

	return 0;
}

