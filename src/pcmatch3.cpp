/*
 * cloudmatch.cpp
 *
 *  Created on: Oct 20, 2019
 *      Author: rob
 */

/**
 * 1. Put each file in a tree, give each a weight.
 * 2. Define a grid.
 * 3. For each cell in the grid, find the points in the neighbourhood for each file,
 *    get the mean of points per file, and the overall weighted mean.
 * 4. Set the weighted mean to the main grid. Set the local mean to the local grid.
 * 5. Smooth the main grid.
 * 6. For each file, find the deviation between the local mean and the smoothed grid.
 * 7. For each point, find the interpolated adjustment and adjust.
 * 8. Write the adjusted files.
 */

#include <string>
#include <vector>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <map>

#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/io/LasReader.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/io/LasWriter.hpp>
#include <pdal/io/LasHeader.hpp>
#include <pdal/Options.hpp>

#include <gdal_priv.h>

#include "ds/mqtree.hpp"
#include "util.hpp"

constexpr double MIN_FLOAT = std::numeric_limits<double>::lowest();
constexpr double MAX_FLOAT = std::numeric_limits<double>::max();

class Point {
private:
	double _x, _y, _z, _w;

public:

	Point() : Point(0, 0, 0) {}

	Point(double x, double y, double z, double weight = 1) :
		_x(x), _y(y), _z(z), _w(weight) {}

	double x() const {
		return _x;
	}

	double y() const {
		return _y;
	}

	double z() const {
		return _z;
	}

	void x(double x) {
		_x = x;
	}

	void y(double y) {
		_y = y;
	}

	void z(double z) {
		_z = z;
	}

	double weight() const {
		return _w;
	}

	void weight(double w) {
		_w = w;
	}

};

class PointFile {
public:
	geo::ds::mqtree<Point> tree;	// Tree to store points from this file.
	std::vector<double> grid;		// The grid of local means for this file.
	std::string file;				// The point source file.
	std::string projection;			// The projection of the las file.
	double minx, miny, maxx, maxy; 	// Boundaries.
	double res;						// The grid size.
	double weight;					// The weight given to the points from this file in the mean (0-1).
	int cols, rows;					// The grid size.

	PointFile(const std::string& file, double res, double weight) :
		file(file),
		minx(MAX_FLOAT), miny(MAX_FLOAT), maxx(MIN_FLOAT), maxy(MIN_FLOAT),
		res(res), weight(weight),
		cols(0), rows(0) {}

	void init() {
		pdal::Option opt("filename", file);
		pdal::Options opts;
		opts.add(opt);
		pdal::PointTable table;
		pdal::LasReader rdr;
		rdr.setOptions(opts);
		rdr.prepare(table);
		pdal::PointViewSet viewSet = rdr.execute(table);
		pdal::PointViewPtr view = *viewSet.begin();
		pdal::Dimension::IdList dims = view->dims();
		pdal::LasHeader hdr = rdr.header();
		projection = hdr.srs().getWKT();

		using namespace pdal::Dimension;

		for (pdal::PointId idx = 0; idx < view->size(); ++idx) {
			double x = view->getFieldAs<double>(Id::X, idx);
			double y = view->getFieldAs<double>(Id::Y, idx);
			if(x < minx) minx = x;
			if(y < miny) miny = y;
			if(x > maxx) maxx = x;
			if(y > maxy) maxy = y;
		}

		minx = std::floor(minx / res) * res;
		miny = std::floor(miny / res) * res;
		maxx = std::floor(maxx / res) * res + res;
		maxy = std::floor(maxy / res) * res + res;

		cols = (int) std::ceil((maxx - minx) / res);
		rows = (int) std::ceil((maxy - miny) / res);

		grid.resize(cols * rows);
		std::fill(grid.begin(), grid.end(), -9999.0);

		tree.init(minx, miny, maxx, maxy);
	}

	void buildTree() {
		std::cout << "Building tree for: " << file << "\n";
		pdal::Option opt("filename", file);
		pdal::Options opts;
		opts.add(opt);
		pdal::PointTable table;
		pdal::LasReader rdr;
		rdr.setOptions(opts);
		rdr.prepare(table);
		pdal::PointViewSet viewSet = rdr.execute(table);
		pdal::PointViewPtr view = *viewSet.begin();
		pdal::Dimension::IdList dims = view->dims();
		pdal::LasHeader hdr = rdr.header();

		using namespace pdal::Dimension;

		double x, y, z;
		int c;

		for (pdal::PointId idx = 0; idx < view->size(); ++idx) {
			c = view->getFieldAs<int>(Id::Classification, idx);
			if(c == 2) {
				x = view->getFieldAs<double>(Id::X, idx);
				y = view->getFieldAs<double>(Id::Y, idx);
				z = view->getFieldAs<double>(Id::Z, idx);
				tree.add(Point(x, y, z));
			}
		}
	}

	void clearTree() {
		tree.clear();
	}

	void set(double x, double y, double v) {
		int col = (int) (x - minx) / res;
		int row = (int) (y - miny) / res;
		if(!(col < 0 || row < 0 || col >= cols || row >= rows))
			grid[row * cols + col] = v;
	}

	double get(double x, double y) {
		int col = (int) (x - minx) / res;
		int row = (int) (y - miny) / res;
		if(!(col < 0 || row < 0 || col >= cols || row >= rows)) {
			return grid[row * cols + col];
		} else  {
			return -9999.0;
		}
	}

};

void usage() {
	std::cout << "Usage: pcmatch3 <pointfile> <outfile> <lasfile [lasfile [...]]>\n";
}

void saveGrid(const std::string& file, const std::vector<double> grid, int cols, int rows, double minx, double miny, double res, const std::string& proj) {
	GDALAllRegister();
	GDALDriverManager* dm = GetGDALDriverManager();
	GDALDriver* drv = dm->GetDriverByName("GTiff");
	GDALDataset* ds = drv->Create(file.c_str(), cols, rows, 1, GDT_Float32, 0);
	double trans[] = {minx, res, 0, miny, 0, res};
	ds->SetGeoTransform(trans);
	ds->SetProjection(proj.c_str());
	if(CE_None != ds->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, cols, rows, (void*) grid.data(), cols, rows, GDT_Float64, 0, 0, 0))
		std::cerr << "Failed to write to raster.\n";
	ds->GetRasterBand(1)->SetNoDataValue(-9999.0);
	GDALClose(ds);
}

double bary(const Point& p, const Point& a, const Point& b, const Point& c) {
	double w1, w2, w3;
	if(a.z() == -9999.0) {
		w1 = 0;
	} else {
		w1 = ((b.y() - c.y()) * (p.x() - c.x()) + (c.x() - b.x()) * (p.y() - c.y())) /
			((b.y() - c.y()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.y() - c.y()));
	}
	if(b.z() == -9999.0) {
		w2 = 0;
	} else {
		w2 = ((c.y() - a.y()) * (p.x() - c.x()) + (a.x() - c.x()) * (p.y() - c.y())) /
			((b.y() - c.y()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.y() - c.y()));
	}
	if(c.z() == -9999.0) {
		w3 = 0;
	} else {
		w3 = 1 - w2 - w1;
	}
	if(w1 + w2 + w3 == 0) {
		return std::nan("");
	} else {
		return a.z() * w1 + b.z() * w2 + c.z() * w3;
	}
}

void gauss2d(std::vector<double>& arr, int size, double sd) {
	arr.resize(std::pow(size * 2 + 1, 2));
	//double a = 1.0 / (2.0 * M_PI * sd * sd);
	double sum = 0;
	for(int r = -size; r < size + 1; ++r) {
		for(int c = -size; c < size + 1; ++c)
			sum += (arr[(r + size) * (size * 2 + 1) + (c + size)] = std::exp(-double(c * c + r * r) / (2.0 * sd * sd)));
	}
	if(sum != 1) {
		for(double& v : arr)
			v /= sum;
	}
}

void convolve(const std::vector<double>& data, int cols, int rows,
		const std::vector<double>& kernel, int size,
		std::vector<double>& smoothed,
		int col = -1, int row = -1) {

	if(col == -1 && row == -1) {
		double v, w, sum, wt;
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				sum = 0;
				wt = 0;
				for(int rr = -size; rr < size + 1; ++rr) {
					for(int cc = -size; cc < size + 1; ++cc) {
						int ccc = cc + c;
						int rrr = rr + r;
						if(ccc < 0 || rrr < 0 || ccc >= cols || rrr >= rows) continue;
						if((v = data[rrr * cols + ccc]) != -9999.0) {
							w = kernel[(rr + size) * (size * 2 + 1) + (cc + size)];
							sum += v * w;
							wt += w;
						}
					}
				}
				if(wt > 0)
					smoothed[r * cols + c] = sum / wt;	// Divide by wt to help with edge effects. Usually it's 1.
			}
		}
	} else {
		double v, w, sum = 0, wt = 0;
		for(int r = -size; r < size + 1; ++r) {
			for(int c = -size; c < size + 1; ++c) {
				int cc = col + c;
				int rr = row + r;
				if(cc < 0 || rr < 0 || cc >= cols || rr >= rows) continue;
				if((v = data[rr * cols + cc]) != -9999.0) {
					w = kernel[(r + size) * (size * 2 + 1) + (c + size)];
					sum += v * w;
					wt += w;
				}
			}
		}
		if(wt > 0)
			smoothed[row * cols + col] = sum / wt;	// Divide by wt to help with edge effects.
	}
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	std::cout << "Creating input entities.\n";
	std::string ptfile = argv[1];
	std::string outfile = argv[2];
	std::vector<PointFile> infiles;
	for(int i = 3; i < argc; ++i)
		infiles.emplace_back(argv[i], 20, 1);

	std::cout << "Building trees.\n";
	for(PointFile& pf : infiles) {
		pf.init();
		pf.buildTree();
	}

	std::vector<std::pair<double, double>> pts;
	{
		std::ifstream ptstr(ptfile);
		std::string line, part;
		std::stringstream ss;
		std::getline(ptstr, line);	// Header
		while(std::getline(ptstr, line)) {
			ss << line;
			std::getline(ss, part, ',');
			double x = atof(part.c_str());
			std::getline(ss, part, ',');
			double y = atof(part.c_str());
			pts.emplace_back(x, y);
			ss.str("");
			ss.clear();
		}
	}

	std::ofstream output(outfile, std::ios::app);
	output << std::setprecision(4) << std::fixed << "file,x,y,z\n";

	double radius = 10;
	std::vector<Point> tres;
	auto iter = std::back_inserter(tres);
	for(const auto& pt : pts) {
		for(PointFile& pf : infiles) {
			pf.tree.search(Point(pt.first, pt.second, 0), radius, iter);
			double s = 0, w = 0;
			for(const Point& rp : tres) {
				double d = std::pow(rp.x() - pt.first, 2.0) + std::pow(rp.y() - pt.second, 2.0);
				if(d == 0) {
					s = rp.z();
					w = 1;
					break;
				} else {
					d = 1 / d;
					s += rp.z() * d;
					w += d;
				}
			}
			tres.clear();
			if(w > 0)
				output << pf.file << "," << pt.first << "," << pt.second << "," << (s/w) << "\n";
		}
	}

	std::cout << "Done.\n";
}
