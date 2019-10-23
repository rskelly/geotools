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
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/io/LasReader.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/io/LasWriter.hpp>
#include <pdal/io/LasHeader.hpp>
#include <pdal/Options.hpp>

#include "ds/mqtree.hpp"
#include "util.hpp"

constexpr double MIN_FLOAT = std::numeric_limits<double>::lowest();
constexpr double MAX_FLOAT = std::numeric_limits<double>::max();

class Point {
public:
	double _x, _y, _z;

	Point() : Point(0, 0, 0) {}

	Point(double x, double y, double z) :
		_x(x), _y(y), _z(z) {}

	double x() const {
		return _x;
	}

	double y() const {
		return _y;
	}

	double z() const {
		return _z;
	}

};

class PointFile {
public:
	std::string file;				// The point source file.
	double minx, miny, maxx, maxy; 	// Boundaries.
	double res;						// The grid size.
	std::vector<double> grid;		// The grid of local means for this file.
	int cols, rows;					// The grid size.
	geo::ds::mqtree<Point> tree;
	double weight;
	std::string projection;

	PointFile(const std::string& file, double res, double weight) :
		file(file),
		minx(MAX_FLOAT), miny(MAX_FLOAT), maxx(MIN_FLOAT), maxy(MIN_FLOAT), res(res),
		cols(0), rows(0),
		weight(weight) {}

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
	std::cout << "Usage: cloudmatch <outfile> <infile [infile [...]]>\n";
}

int gridNo = 0;
#include <gdal_priv.h>

void saveGrid(const std::vector<double> grid, int cols, int rows, double minx, double miny, double res, const std::string& proj) {
	std::string file = "grid_" + std::to_string(++gridNo) + ".tif";
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

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	double res = 100;

	std::cout << "Creating input entities.\n";
	std::string outdir = argv[1];
	std::vector<PointFile> infiles;
	for(int i = 2; i < argc; ++i)
		infiles.emplace_back(argv[i], res, 1);

	std::cout << "Initializing entities and bounds.\n";
	double minx = MAX_FLOAT;
	double miny = MAX_FLOAT;
	double maxx = MIN_FLOAT;
	double maxy = MIN_FLOAT;
	for(PointFile& pf : infiles) {
		pf.init();
		if(minx > pf.minx) minx = pf.minx;
		if(miny > pf.miny) miny = pf.miny;
		if(maxx < pf.maxx) maxx = pf.maxx;
		if(maxy < pf.maxy) maxy = pf.maxy;
	}

	minx = std::floor(minx / res) * res;
	miny = std::floor(miny / res) * res;
	maxx = std::floor(maxx / res) * res + res;
	maxy = std::floor(maxy / res) * res + res;

	std::cout << "Building trees.\n";
	for(PointFile& pf : infiles)
		pf.buildTree();

	int cols = (int) std::ceil((maxx - minx) / res);
	int rows = (int) std::ceil((maxy - miny) / res);

	std::vector<double> grid(cols * rows);
	std::fill(grid.begin(), grid.end(), -9999.0);

	Point pt(0, 0, 0);
	std::vector<Point> pts;
	auto ins = std::back_inserter(pts);

	std::cout << "Creating grids.\n";
	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {
			double x = minx + col * res + res * 0.5;
			double y = miny + row * res + res * 0.5;

			pt._x = x;
			pt._y = y;

			double sum = 0;
			double w = 0;

			for(PointFile& pf : infiles) {
				size_t count = pf.tree.search(pt, res, ins);
				if(count) {
					double s = 0;
					for(const Point& p : pts)
						s += p.z();
					s /= count;
					sum += s * pf.weight;
					w += pf.weight;
					pf.set(x, y, s);
				}
				pts.clear();
			}

			if(w)
				grid[row * cols + col] = sum / w;
		}
	}

	std::cout << "Saving grids.\n";
	for(PointFile& pf : infiles)
		saveGrid(pf.grid, pf.cols, pf.rows, pf.minx, pf.miny, res, pf.projection);
	saveGrid(grid, cols, rows, minx, miny, res, infiles.front().projection);

	std::cout << "Clearing trees.\n";
	for(PointFile& pf : infiles)
		pf.clearTree();

	// now smooth the main grid
	std::cout << "Smoothing the grid.\n";
	{
		int ksize = 11;
		std::vector<std::tuple<int, int, double>> kernel(ksize * ksize);
		for(int r = -ksize / 2; r < ksize / 2 + 1; ++r) {
			for(int c = -ksize / 2; c < ksize / 2 + 1; ++c)
				kernel[(r + ksize / 2) * ksize + (c + ksize / 2)] = std::make_tuple(c, r, 1.0 - std::min(1.0, (c * c + r * r) / std::pow(ksize / 2, 2)));
		}
		std::vector<double> smoothed(cols * rows);
		std::fill(smoothed.begin(), smoothed.end(), -9999.0);
		for(int row = 0; row < rows; ++row) {
			for(int col = 0; col < cols; ++col) {
				double v, sum = 0, weight = 0;
				int ct = 0;
				for(const auto& it : kernel) {
					int cc = col + std::get<0>(it);
					int rr = row + std::get<1>(it);
					double w = std::get<2>(it);
					if(!(cc < 0 || rr < 0 || cc >= cols || rr >= rows)
							&& (v = grid[rr * cols + cc]) != -9999.0) {
						sum += v * w;
						weight += w;
						++ct;
					}
				}

				if(ct)
					smoothed[row * cols + col] = sum / weight;
			}
		}
		grid.swap(smoothed);
	}

	// now adjust and write each file grid using the deviation.
	std::cout << "Adjusting file grids.\n";
	for(PointFile& pf : infiles) {
		for(int pr = 0; pr < pf.rows; ++pr) {
			for(int pc = 0; pc < pf.cols; ++pc) {
				double pv = pf.grid[pr * pf.cols + pc];
				if(pv != -9999.0) {
					double px = pf.minx + pc * res + res * 0.5;
					double py = pf.miny + pr * res + res * 0.5;
					int gc = (int) (px - minx) / res;
					int gr = (int) (py - miny) / res;
					double gv = grid[gr * cols + gc];
					if(gv == -9999.0)
						throw std::runtime_error("Found illegal nodata in grid.");
					pf.grid[pr * pf.cols + pc] = gv - pv;
				}
			}
		}
	}

	std::cout << "Saving grids.\n";
	for(PointFile& pf : infiles)
		saveGrid(pf.grid, pf.cols, pf.rows, pf.minx, pf.miny, res, pf.projection);
	saveGrid(grid, cols, rows, minx, miny, res, infiles.front().projection);

	// now adjust each las file and write new ones.
	std::cout << "Adjusting files...\n";

	geo::util::makedir(outdir);

	for(PointFile& pf : infiles) {
		pdal::Option opt("filename", pf.file);
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

		double minz = MAX_FLOAT;
		double maxz = MIN_FLOAT;

		double x, y, z, pv;
		int col, row, col0, row0;
		double cx, cy, cx0, cy0, cx1, cy1;

		for (pdal::PointId idx = 0; idx < view->size(); ++idx) {
			x = view->getFieldAs<double>(Id::X, idx);
			y = view->getFieldAs<double>(Id::Y, idx);
			z = view->getFieldAs<double>(Id::Z, idx);

			col = (int) (x - pf.minx) / res;
			row = (int) (y - pf.miny) / res;

			// The centre of the cell containing the point.
			cx = pf.minx + col * res + res * 0.5;
			cy = pf.miny + row * res + res * 0.5;

			col0 = col;
			if(x < cx) {
				col0 = col > 0 ? col - 1 : col;
			} else if(x > cx) {
				col0 = col < cols - 1 ? col + 1 : col;
			}

			row0 = row;
			if(y < cy) {
				row0 = row > 0 ? row - 1 : row;
			} else if(y > cy) {
				row0 = row < rows - 1 ? row + 1 : row;
			}

			cx0 = pf.minx + col0 * res + res * 0.5;
			cy0 = pf.miny + row * res + res * 0.5;

			cx1 = pf.minx + col * res + res * 0.5;
			cy1 = pf.miny + row0 * res + res * 0.5;

			pv = bary(Point(x, y, 0),
				Point(cx, cy, pf.grid[row * pf.cols + col]),
				Point(cx0, cy0, pf.grid[row * pf.cols + col0]),
				Point(cx1, cy1, pf.grid[row0 * pf.cols + col])
			);

			if(std::isnan(pv))
				pv = pf.grid[row * pf.cols + col];

			z += pv;

			view->setField(Id::Z, idx, z);

			if(z < minz) minz = z;
			if(z > maxz) maxz = z;
		}

		std::string outfile = geo::util::join(outdir, geo::util::basename(pf.file) + "_adj.las");
		if(geo::util::isfile(outfile))
			geo::util::rem(outfile);
		opt = pdal::Option("filename", outfile);
		opts.replace(opt);
		pdal::LasWriter wtr;
		pdal::BufferReader brdr;
		brdr.addView(view);
		wtr.setInput(brdr);
		wtr.setOptions(opts);
		wtr.setSpatialReference(rdr.getSpatialReference());

		wtr.prepare(table);
		wtr.execute(table);
	}

	std::cout << "Done.\n";
}
