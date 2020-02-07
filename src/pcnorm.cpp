/*
 * pcnorm2.cpp
 *
 *  Created on: May 5, 2019
 *      Author: rob
 */

#include <limits>
#include <fstream>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <fstream>

#include "pointcloud_reader.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3   Point;
typedef Gt::Point_2  Point2;
typedef CGAL::Delaunay_triangulation_2<Gt>::Face_handle FaceHandle;

size_t buildGrid(const std::string& infile, double* bounds, double res, int& cols, int& rows, std::vector<double>& grid) {

	for(int i = 0; i < 3; ++i) bounds[i] = G_DBL_MAX_POS;
	for(int i = 3; i < 6; ++i) bounds[i] = G_DBL_MAX_NEG;

	pointcloud_reader rdr;

	std::cout << "Computing bounds\n";
	rdr.open(infile, true);
	rdr.extendBounds(bounds);

	// Get grid size.
	cols = (int) std::ceil((bounds[3] - bounds[0]) / res);
	rows = (int) std::ceil((bounds[4] - bounds[1]) / res);

	// Accumulate the cell elevations from points.
	grid.resize(cols * rows);
	std::fill(grid.begin(), grid.end(), 0);
	std::vector<int> counts(cols * rows);
	std::fill(counts.begin(), counts.end(), 0);
	double x, y, z;
	int cls;

	while(rdr.next(x, y, z, cls)) {
		if(cls == 2) {
			int c = (int) (x - bounds[0]) / res;
			int r = (int) (y - bounds[1]) / res;
			if(c < 0 || r < 0 || c >= cols || r >= rows)
				continue; //throw std::runtime_error("out of range");
			grid[r * cols + c] += z;
			counts[r * cols + c]++;
		}
	}
	rdr.close();

	size_t count = 0;

	// Compute the means.
	for(size_t i = 0; i < grid.size(); ++i) {
		if(counts[i]) {
			grid[i] /= counts[i];
			++count;
		} else {
			grid[i] = -9999;
		}
	}

	size_t filled;
	double s, v;
	int ct;
	do {
		filled = 0;
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				s = 0;
				ct = 0;
				if(grid[r * cols + c] == -9999) {
					for(int rr = r - 1; rr < r + 2; ++rr) {
						for(int cc = c - 1; cc < c + 2; ++cc) {
							if(!(rr < 0 || cc < 0 || (rr == r && cc == c) || rr >= rows || cc >= cols)
									&& (v = grid[rr * cols + cc]) != -9999) {
								s += v;
								++ct;
							}
						}
					}
				}
				if(ct) {
					grid[r * cols + c] = s / ct;
					++filled;
				}
			}
		}
	} while(filled);

	return count;
}

double bary(double x, double y,
		double x0, double y0, double z0,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	double w0 = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) / ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2));
	double w1 = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) / ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2));
	double w2 = 1 - w0 - w1;
	return (w0 * z0) + (w1 * z1) + (w2 * z2);
}

bool buildMesh(std::vector<double>& grid, int cols, int rows, double* bounds, double res, Delaunay& mesh) {

	std::vector<Point> pts;
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {
			double x = bounds[0] + c * res + res * 0.5;
			double y = bounds[1] + r * res + res * 0.5;
			pts.emplace_back(x, y, grid[r * cols + c]);
		}
	}
	mesh = Delaunay(pts.begin(), pts.end());
	std::cout << "Delaunay: " << mesh.number_of_vertices() << "vertices.\n";
	return true;
}

bool normalize(const std::string& infile, Delaunay& mesh, const std::string& outfile) {

	pointcloud_reader rdr;

	double x, y, z;
	int cls;
	FaceHandle faceh;

	rdr.open(infile, false);

	std::cout << std::setprecision(3) << std::fixed;
	size_t n = 0;
	while(rdr.next(x, y, z, cls)) {
		Point2 pt(x, y, z);
		faceh = mesh.locate(pt, faceh);
		auto v1 = faceh->vertex(0)->point();
		auto v2 = faceh->vertex(1)->point();
		auto v3 = faceh->vertex(2)->point();
		double zi = bary(x, y, v1.x(), v1.y(), v1.z(), v2.x(), v2.y(), v2.z(), v3.x(), v3.y(), v3.z());
		rdr.setZ(zi);
		if(++n % 1000 == 0)
			std::cout << n << " of " << rdr.size() << " (" << ((float) n / rdr.size()) << "%)\n";
	}

	rdr.save(outfile);

	return true;
}

bool normalizeRaster(const std::string& infile, std::vector<double>& grid, double* bounds,
		double res, int cols, int rows, const std::string& outfile) {

	pointcloud_reader rdr;

	double x, y, z;
	int cls;

	rdr.open(infile, false);

	size_t n = 0;
	size_t idx;
	while(rdr.next(x, y, z, cls, idx)) {
		int col = (int) (x - bounds[0]) / res;
		int row = (int) (y - bounds[1]) / res;
		int col0 = col;
		int row0 = row;
		double cx0 = bounds[0] + col * res + 0.5 * res;
		double cy0 = bounds[1] + row * res + 0.5 * res;
		double z0 = grid[row * cols + col];
		double z1 = -9999;
		double z2 = -9999;
		int ci = x < cx0 ? -1 : 1;
		int cl = x < cx0 ? 0 : cols - 1;
		int ri = y < cy0 ? -1 : 1;
		int rl = y < cy0 ? 0 : rows - 1;
		col0 += ci;
		do {
			if(col0 < 0 || col0 >= cols || (z1 = grid[row * cols + col0]) != -9999)
				break;
			col0 += ci;
		} while(col0 != cl);
		row0 += ri;
		do {
			if(row0 < 0 || row0 >= rows || (z2 = grid[row0 * cols + col]) != -9999)
				break;
			row0 += ri;
		} while(row0 != rl);
		if(z1 == -9999) z1 = z0;
		if(z2 == -9999) z2 = z0;
		double cx1 = bounds[0] + col0 * res + 0.5 * res;
		double cy1 = bounds[1] + row * res + 0.5 * res;
		double cx2 = bounds[0] + col * res + 0.5 * res;
		double cy2 = bounds[1] + row0 * res + 0.5 * res;
		double zi = bary(x, y, cx0, cy0, z0, cx1, cy1, z1, cx2, cy2, z2);
		rdr.setZ(idx, z - zi);
		if(++n % 1000000 == 0)
			std::cout << (n * 100) / rdr.size() << "%\n";
	}

	rdr.save(outfile);

	return true;
}

int main(int argc, char** argv) {

	if(argc < 4) {
		std::cerr << "Usage: pcnorm <outfile (.las)> <resolution> <infile(s) (.las)>\n";
		return 1;
	}

	std::string outfile = argv[1];
	double resolution = atof(argv[2]);
	std::string infile = argv[3];

	double bounds[6];
	int cols, rows;
	std::vector<double> grid;
	Delaunay mesh;

	std::cout << "Building grid.\n";
	if(!buildGrid(infile, bounds, resolution, cols, rows, grid)) {
		std::cerr << "Failed to build grid.\n";
		return 1;
	}

	//std::cout << "Building mesh.\n";
	//if(!buildMesh(grid, cols, rows, bounds, resolution, mesh)) {
	//	std::cerr << "Failed to build mesh.\n";
	//	return 1;
	//}
	//grid.clear();

	std::cout << "Normalizing.\n";
	//if(!normalize(infile, mesh, outfile)) {
	if(!normalizeRaster(infile, grid, bounds, resolution, cols, rows, outfile)) {
		std::cerr << "Failed to normalize.\n";
		return 1;
	}

	return 0;
}
