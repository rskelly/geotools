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
	double _x, _y, _z;

public:

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

	void x(double x) {
		_x = x;
	}

	void y(double y) {
		_y = y;
	}

	void z(double z) {
		_z = z;
	}

};

class PointFile {
public:
	std::string file;				// The point source file.
	double minx, miny, maxx, maxy; 	// Boundaries.
	double res;						// The grid size.
	std::vector<double> grid;		// The grid of local means for this file.
	int cols, rows;					// The grid size.
	geo::ds::mqtree<Point> tree;	// Tree to store points from this file.
	double weight;					// The weight given to the points from this file in the mean (0-1).
	std::string projection;			// The projection of the las file.

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

void gauss(std::vector<double>& arr, int size, double sd) {
	arr.resize(std::pow(size * 2 + 1, 2));
	double a = 1.0 / (2.0 * M_PI * sd * sd);
	for(int r = -size; r < size + 1; ++r) {
		for(int c = -size; c < size + 1; ++c) {
			arr[(r + size) * (size * 2 + 1) + (c + size)] = a * std::exp(-double(c * c + r * r) / (2.0 * sd * sd));
		}
	}
}

void convolve(std::vector<double>& data, int cols, int rows, const std::vector<double>& kernel, int size, int col = -1, int row = -1) {
	if(col == -1 && row == -1) {
		double v, w, sum = 0, wt = 0;
		std::vector<double> smoothed(data.size());
		std::fill(smoothed.begin(), smoothed.end(), -9999.0);
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				for(int rr = r - size; rr < r + size + 1; ++rr) {
					for(int cc = c - size; cc < c + size + 1; ++cc) {
						if(cc < 0 || rr < 0 || cc >= cols || rr >= rows) continue;
						if((v = data[rr * cols + cc]) != -9999.0) {
							w = kernel[(rr + size) * (size * 2 + 1) + (cc + size)];
							sum += v * w;
							wt += w;
						}
					}
				}
				smoothed[r * cols + c] = sum / wt;	// Divide by wt to help with edge effects.
			}
		}
		data.swap(smoothed);
	} else {
		double v, w, sum = 0, wt = 0;
		for(int rr = row - size; rr < row + size + 1; ++rr) {
			for(int cc = col - size; cc < col + size + 1; ++cc) {
				if(cc < 0 || rr < 0 || cc >= cols || rr >= rows) continue;
				if((v = data[rr * cols + cc]) != -9999.0) {
					w = kernel[(rr + size) * (size * 2 + 1) + (cc + size)];
					sum += v * w;
					wt += w;
				}
			}
		}
		data[row * cols + col] = sum / wt;	// Divide by wt to help with edge effects.
	}
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	double res = 20;

	for(int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-r") {
			res = atof(argv[++i]);
		}
	}

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

			pt.x(x);
			pt.y(y);

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

	//std::cout << "Saving grids.\n";
	//for(PointFile& pf : infiles)
	//	saveGrid(pf.grid, pf.cols, pf.rows, pf.minx, pf.miny, res, pf.projection);
	//saveGrid(grid, cols, rows, minx, miny, res, infiles.front().projection);

	std::cout << "Clearing trees.\n";
	for(PointFile& pf : infiles)
		pf.clearTree();


	// now smooth the main grid
	std::cout << "Smoothing the grid.\n";
	{
		int crad = 10;	// Radius of kernel in cells (total width, 21).

		// Collect the variances within the kernel region. These are used to approximate the slope-ness
		// of the terrain. TODO: Use an actual slope.
		std::cout << "Collecting stats.\n";
		std::vector<double> quantiles(crad); 	// Stores the boundaries for the equal-count partition of unique variance values.
		std::vector<double> stats(cols * rows);
		std::fill(stats.begin(), stats.end(), -9999.0);
		{
			std::set<double> variances;				// Stores the unique values of the variances.
			for(int row = 0; row < rows; ++row) {
				for(int col = 0; col < cols; ++col) {
					double v, mean = 0;
					int ct = 0;
					for(int rr = row - crad; rr < row + crad + 1; ++rr) {
						for(int cc = col - crad; cc < col + crad + 1; ++cc) {
							if(cc < 0 || rr < 0 || cc >= cols || rr >= rows) continue;
							if((v = grid[rr * cols + cc]) != -9999.0) {
								mean += v;
								++ct;
							}
						}
					}

					if(!ct) continue;

					mean /= ct;

					double var = 0;
					for(int rr = row - crad; rr < row + crad + 1; ++rr) {
						for(int cc = col - crad; cc < col + crad + 1; ++cc) {
							if(cc < 0 || rr < 0 || cc >= cols || rr >= rows) continue;
							if((v = grid[rr * cols + cc]) != -9999.0)
								var += std::pow(v - mean, 2.0);
						}
					}

					var /= ct;

					stats[row * cols + col] = var;

					variances.insert(var);
				}
			}

			// Get the equal-area boundaries of the variance set.
			int step = variances.size() / crad;
			std::vector<double> variances0(variances.begin(), variances.end());
			for(int i = 0; i < crad; ++i)
				quantiles[i] = variances0[(i + 1) * step];

			// Smooth the stats raster. There's a lot of noise here.
			{
				std::vector<double> gf;
				gauss(gf, 1, 0.5);
				convolve(stats, cols, rows, gf, 1);
			}

			saveGrid("stats.tif", stats, cols, rows, minx, miny, res, infiles.front().projection);
		}

		// Make the smoothing kernels. There are {crad} kernels, one for each
		// radius from 1 to {crad} inclusive. The kernels are distance-weighted averages.
		std::vector<std::vector<double>> kernels(crad);
		{
			std::vector<double> gf;
			for(int i = 1; i <= crad; ++i)
				gauss(kernels[i - 1], i, i * 0.5);
		}

		// Smooth the raster.
		std::vector<double> smoothed(cols * rows);
		std::fill(smoothed.begin(), smoothed.end(), -9999.0);
		for(int row = 0; row < rows; ++row) {
			for(int col = 0; col < cols; ++col) {

				// Get the correct kernel for the variance for this cell.
				// The lowest variance gets the largest kernel.
				double var = stats[row * cols + col];
				if(var == -9999.0) continue;
				int varidx = 0;
				for(double q : quantiles) {
					if(q < var)
						++varidx;
				}
				convolve(grid, cols, rows, kernels[varidx], varidx + 1, col, row);
			}
		}

		saveGrid("grid.tif", grid, cols, rows, minx, miny, res, infiles.front().projection);
	}

	if(false) {
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

		//std::cout << "Saving grids.\n";
		//for(PointFile& pf : infiles)
		//	saveGrid(pf.grid, pf.cols, pf.rows, pf.minx, pf.miny, res, pf.projection);
		//saveGrid(grid, cols, rows, minx, miny, res, infiles.front().projection);

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

			double x, y, z, cx, cy;
			double d, v, w, ss, sw;
			bool skip;
			int col, row;

			for (pdal::PointId idx = 0; idx < view->size(); ++idx) {
				x = view->getFieldAs<double>(Id::X, idx);
				y = view->getFieldAs<double>(Id::Y, idx);
				z = view->getFieldAs<double>(Id::Z, idx);

				col = (int) (x - pf.minx) / res;
				row = (int) (y - pf.miny) / res;

				ss = 0;
				sw = 0;
				skip = false;
				for(int rr = row - 1; !skip && rr < row + 2; ++rr) {
					for(int cc = col - 1; !skip && cc < col + 2; ++cc) {
						if(cc < 0 || rr < 0 || rr >= pf.rows || cc >= pf.cols) continue;
						if((v = pf.grid[rr * pf.cols + cc]) != -9999.0) {
							cx = pf.minx + cc * res + res * 0.5;
							cy = pf.miny + rr * res + res * 0.5;
							d = std::pow(cx - x, 2.0) + std::pow(cy - y, 2.0);
							if(d == 0) {
								ss = v;
								sw = 1.0;
								skip = true;
								break;
							} else {
								w = 1.0 - std::min(1.0, d / (res * res));
								ss += v * w;
								sw += w;
							}
						}
					}
				}

				if(sw != 0) {
					z += ss / sw;
				} else {
					std::cout << ss << ", " << sw << ", " << col << ", " << row << "\n";
				}

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
	}
	std::cout << "Done.\n";
}
