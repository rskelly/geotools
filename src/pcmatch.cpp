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

#include <liblas/liblas.hpp>

#include "ds/mqtree.hpp"
#include "util.hpp"

constexpr float MIN_FLOAT = std::numeric_limits<float>::lowest();
constexpr float MAX_FLOAT = std::numeric_limits<float>::max();

class Point {
public:
	float _x, _y, _z;

	Point() : Point(0, 0, 0) {}

	Point(float x, float y, float z) :
		_x(x), _y(y), _z(z) {}

	float x() const {
		return _x;
	}

	float y() const {
		return _y;
	}

	float z() const {
		return _z;
	}

};

class PointFile {
public:
	std::string file;			// The point source file.
	float minx, miny, maxx, maxy; 			//
	float res;					// The grid size.
	std::vector<float> grid;	// The grid of local means for this file.
	int cols, rows;				// The grid size.
	geo::ds::mqtree<Point> tree;
	float weight;

	PointFile(const std::string& file, float res, float weight) :
		file(file),
		minx(MAX_FLOAT), miny(MAX_FLOAT), maxx(MIN_FLOAT), maxy(MIN_FLOAT), res(res),
		cols(0), rows(0),
		weight(weight) {}

	void init() {
		liblas::ReaderFactory rf;
		std::ifstream str(file);
		liblas::Reader rdr = rf.CreateWithStream(str);
		liblas::Header hdr = rdr.GetHeader();
		while(rdr.ReadNextPoint()) {
			const liblas::Point& pt = rdr.GetPoint();
			float x = (float) pt.GetX();
			float y = (float) pt.GetY();
			if(x < minx) minx = x;
			if(y < miny) miny = y;
			if(x > maxx) maxx = x;
			if(y > maxy) maxy = y;
		}
		cols = (int) std::ceil((maxx - minx) / res);
		rows = (int) std::ceil((maxy - miny) / res);
		grid.resize(cols * rows);
		tree.init(minx, miny, maxx, maxy);
	}

	void buildTree() {
		liblas::ReaderFactory rf;
		std::ifstream str(file);
		liblas::Reader rdr = rf.CreateWithStream(str);
		liblas::Header hdr = rdr.GetHeader();
		while(rdr.ReadNextPoint()) {
			const liblas::Point& pt = rdr.GetPoint();
			if(pt.GetClassification().GetClass() == 2)
				tree.add(Point(pt.GetX(), pt.GetY(), pt.GetZ()));
		}
	}

	void clearTree() {
		tree.clear();
	}

	void set(float x, float y, float v) {
		int col = (int) (x - minx) / res;
		int row = (int) (y - miny) / res;
		grid[row * rows + col] = v;
	}

	float get(float x, float y) {
		int col = (int) (x - minx) / res;
		int row = (int) (y - miny) / res;
		return grid[row * rows + col];
	}

};

void usage() {
	std::cout << "Usage: cloudmatch <outfile> <infile [infile [...]]>\n";
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	float res = 100;

	std::cout << "Creating input entities.\n";
	std::string outdir = argv[1];
	std::vector<PointFile> infiles;
	for(int i = 2; i < argc; ++i)
		infiles.emplace_back(argv[i], res, 1);

	std::cout << "Initializing entities and bounds.\n";
	float minx = std::numeric_limits<float>::max();
	float miny = std::numeric_limits<float>::max();
	float maxx = std::numeric_limits<float>::lowest();
	float maxy = std::numeric_limits<float>::lowest();
	for(PointFile& pf : infiles) {
		pf.init();
		if(minx > pf.minx) minx = pf.minx;
		if(miny > pf.miny) miny = pf.miny;
		if(maxx < pf.maxx) maxx = pf.maxx;
		if(maxy < pf.maxy) maxy = pf.maxy;
	}

	std::cout << "Building trees.\n";
	for(PointFile& pf : infiles)
		pf.buildTree();

	int cols = (int) std::ceil((maxx - minx) / res);
	int rows = (int) std::ceil((maxy - miny) / res);

	std::vector<float> grid(cols * rows);

	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {
			float x = minx + col * res + res * 0.5;
			float y = miny + row * res + res * 0.5;

			Point pt(x, y, 0);
			std::vector<Point> pts;
			auto ins = std::back_inserter(pts);

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
				} else {
					pf.set(x, y, -9999.0);
				}
				pts.clear();
			}

			if(w) {
				grid[row * rows + col] = sum / w;
			} else {
				grid[row * rows + col] = -9999.0;
			}
		}
	}

	std::cout << "Clearing trees.\n";
	for(PointFile& pf : infiles)
		pf.clearTree();

	// now smooth the main grid
	std::cout << "Smoothing the grid.\n";
	for(int row = 0; row < rows; ++row) {
		for(int col = 0; col < cols; ++col) {
			float v, sum = 0;
			int ct = 0;
			for(int rr = row - 5; rr < row + 6; ++rr) {
				for(int cc = col - 5; cc < col + 6; ++cc) {
					if(cc < 0 || rr < 0 || cc >= cols || rr >= rows)
						continue;
					if((v = grid[rr * rows + cc]) != -9999.0) {
						sum += v;
						++ct;
					}
				}
			}
			if(ct) {
				grid[row * rows + col] = sum / ct;
			}
		}
	}

	// now adjust and write each file grid using the deviation.
	std::cout << "Adjusting file grids.\n";
	for(PointFile& pf : infiles) {
		for(int pr = 0; pr < pf.rows; ++pr) {
			for(int pc = 0; pc < pf.cols; ++pc) {
				float pv = pf.grid[pr * pf.rows + pc];
				if(pv != -9999.0) {
					float px = pf.minx + pc * res + res * 0.5;
					float py = pf.miny + pc * res + res * 0.5;
					int gc = (int) (px - minx) / res;
					int gr = (int) (py - miny) / res;
					float gv = grid[gr * rows + gc];
					if(gv != -9999.0)
						pf.grid[pr * pf.rows + pc] = gv - pv;
				}
			}
		}
	}

	// now adjust each las file and write new ones.
	std::cout << "Adjusting files...\n";

	geo::util::makedir(outdir);

	for(PointFile& pf : infiles) {
		liblas::ReaderFactory rf;
		std::ifstream istr(pf.file);
		liblas::Reader rdr = rf.CreateWithStream(istr);
		liblas::Header hdr(rdr.GetHeader());

		std::string outfile = geo::util::join(outdir, geo::util::basename(pf.file) + "_adj.las");
		std::ofstream ostr(outfile);
		liblas::WriterFactory wf;
		liblas::Writer wtr = wf.CreateWithStream(ostr, hdr);

		float minz = MAX_FLOAT;
		float maxz = MIN_FLOAT;
		while(rdr.ReadNextPoint()) {
			const liblas::Point& pt = rdr.GetPoint();
			int col = (int) (pt.GetX() - pf.minx) / res;
			int row = (int) (pt.GetY() - pf.miny) / res;
			float v = pf.grid[row * pf.rows + col];
			if(v != -9999.0) {
				liblas::Point pt0(pt);
				pt0.SetZ(pt.GetZ() + v);
				if(pt0.GetZ() < minz) minz = pt0.GetZ();
				if(pt0.GetZ() > maxz) maxz = pt0.GetZ();
				wtr.WritePoint(pt0);
			}
		}

		hdr.SetMax(pf.maxx, pf.maxy, maxz);
		hdr.SetMin(pf.minx, pf.miny, minz);

		wtr.WriteHeader();
	}

	std::cout << "Done.\n";
}
