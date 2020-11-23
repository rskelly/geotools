/*
 * pcnorm2.cpp
 *
 *  Created on: May 5, 2019
 *      Author: rob
 */

#include <limits>
#include <fstream>

#include "geo.hpp"
#include "util.hpp"
#include "pointcloud.hpp"

using namespace geo::util;

bool buildGrid(const std::vector<std::string>& infiles, double* bounds, double res, int& cols, int& rows, std::vector<float>& grid) {

	// Increase bounds enough to add 2 cells all around.
	bounds[0] -= res;
	bounds[1] += res;
	bounds[2] -= res;
	bounds[3] += res;

	// Get grid size.
	cols = (int) std::ceil((bounds[1] - bounds[0]) / res);
	rows = (int) std::ceil((bounds[3] - bounds[2]) / res);

	std::vector<float> weights(cols * rows);
	size_t count = 0;
	{
		// To compute weighted mean, need weights and values arrays.
		grid.resize(cols * rows);
		std::fill(grid.begin(), grid.end(), 0);
		std::fill(weights.begin(), weights.end(), 0);

		double px, py, pz;
		double x, y, d, w;
		int col, row;
		double rad = std::pow(res, 2);

		size_t num = 0;
		for(const std::string& infile : infiles) {
			std::cout << ++num << " of " << infiles.size() << "\n";
			geo::pc::PCFile rdr(infile);
			geo::pc::Point pt;
			while(rdr.next(pt)) {
				if(pt.classId() != 2)
					continue;

				px = pt.x();
				py = pt.y();
				pz = pt.z();
				col = (int) (px - bounds[0]) / res;
				row = (int) (py - bounds[2]) / res;

				for(int r = row - 1; r < row + 2; ++r) {
					for(int c = col - 1; c < col + 2; ++c) {
						if(c >= 0 && r >= 0 && c < cols && r < rows) {

							x = bounds[0] + (c * res) + res * 0.5;
							y = bounds[2] + (r * res) + res * 0.5;
							d = std::pow(x - px, 2.0) + std::pow(y - py, 2.0);

							if(d > rad)
								continue;

							w = 1.0 - d / rad;

							// Accumulate the weighted heights and weights.
							grid[r * cols + c] += pz * w;
							weights[r * cols + c] += w;
							++count;
						}
					}
				}
			}
		}
	}

	if(!count) {
		std::cout << "There are no ground points. Quitting.\n";
		return false;
	}

	// Normalize by weights.
	for(size_t i = 0; i < grid.size(); ++i) {
		if(weights[i] > 0)
			grid[i] /= weights[i];
	}

	// Fill in zeroes.
	std::cout << "Filling gaps\n";
	for(size_t i = 0; i < grid.size(); ++i) {
		if(i % cols == 0)
			std::cout << "Row " << (i / cols) << " of " << rows << "\n";
		if(weights[i] == 0) {
			double x, y, d, w0, rad, w = 0, s = 0;
			int o = 2;
			int col = i % cols;
			int row = i / cols;
			double cx = bounds[0] + col * res + res * 0.5;
			double cy = bounds[2] + row * res + res * 0.5;
			do {
				rad = std::pow(o * res, 2);
				for(int r = row - o; r < row + o + 1; ++r) {
					for(int c = col - o; c < col + o + 1; ++c) {
						if(c >= 0 && r >= 0 && c < cols && r < rows) {

							if(weights[r * cols + c] == 0)
								continue;

							x = bounds[0] + (c * res) + res * 0.5;
							y = bounds[2] + (r * res) + res * 0.5;
							d = std::pow(x - cx, 2.0) + std::pow(y - cy, 2.0);

							if(d > rad)
								continue;

							w0 = 1.0 - d / rad;

							// Accumulate the weighted heights and weights.
							s += grid[r * cols + c] * w0;
							w += w0;
						}
					}
				}
				++o;
			} while(w == 0);

			grid[i] = s / w;
		}
	}
	return true;
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

void normalize(const std::vector<std::string>& infiles, geo::pc::PCWriter& wtr,
		double* bounds, double res, int cols, int /*rows*/, std::vector<float>& grid) {

	bounds[4] = std::numeric_limits<double>::max();
	bounds[5] = std::numeric_limits<double>::lowest();

	double z, px, py, pz, cx0, cy0, cz0, cx1, cy1, cz1, cx2, cy2, cz2, nz;
	int col, row, col0, row0;

	size_t num = 0;
	for(const std::string& infile : infiles) {
		std::cout << ++num << " of " << infiles.size() << "\n";
		geo::pc::PCFile rdr(infile);
		geo::pc::Point pt;
		while(rdr.next(pt)) {
			px = pt.x();
			py = pt.y();
			pz = pt.z();
			// Point's home cell.
			col = (int) (px - bounds[0]) / res;
			row = (int) (py - bounds[2]) / res;
			// Center of home cell.
			cx0 = bounds[0] + (col * res) + res * 0.5;
			cy0 = bounds[2] + (row * res) + res * 0.5;
			cz0 = grid[row * cols + col];
			// Cell offsets.
			col0 = px < cx0 ? col - 1 : col + 1;
			row0 = py < cy0 ? row - 1 : row + 1;
			// Centers of offset cells.
			cx1 = bounds[0] + (col0 * res) + res * 0.5;
			cy1 = bounds[2] + (row * res) + res * 0.5;
			cz1 = grid[row * cols + col0];
			cx2 = bounds[0] + (col * res) + res * 0.5;
			cy2 = bounds[2] + (row0 * res) + res * 0.5;
			cz2 = grid[row0 * cols + col];

			if(cz0 == -9999.0 || cz1 == -9999.0 || cz2 == -9999.0)
				continue;

			// Get the barycentric z
			nz = bary(px, py, cx0, cy0, cz0, cx1, cy1, cz1, cx2, cy2, cz2);
			// Make point and write it.
			pt.z((z = pz - nz));
			wtr.addPoint(pt);
			// Adjust bounds.
			if(z < bounds[4]) bounds[4] = z;
			if(z > bounds[5]) bounds[5] = z;
		}
	}
}

void usage() {
	std::cerr << "Usage: pcnorm [options] <outfile (.las)> <resolution> <infile(s) (.las)>\n"
			<< " -f   Force overwrite of existing files.\n";

}
int main(int argc, char** argv) {

	if(argc < 4) {
		usage();
		return 1;
	}

	std::string outfile;
	std::vector<std::string> infiles;
	double resolution = 0;
	bool force = false;

	int mode = 0;
	for(int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if(arg == "-f") {
			force = true;
		} else if(mode == 0) {
			++mode;
			outfile = arg;
		} else if(mode == 1) {
			++mode;
			resolution = atof(arg.c_str());
		} else {
			infiles.push_back(arg);
		}
	}

	if(!checkValidInputFiles(infiles)) {
		g_error("At least one of the input files is invalid.");
		usage();
		return 1;
	}

	if(resolution <= 0 || std::isnan(resolution)) {
		g_error("The resolution " << resolution << " is invalid.");
		usage();
		return 1;
	}

	if(outfile.empty()) {
		g_error("Output file not given.");
		usage();
		return 1;
	}

	if(!geo::util::safeToWrite(outfile, force)) {
		g_error("The file " << outfile << " is not safe to write to. "
				<< "If it is a file, use the -f flag. If it is a directory, "
				<< "choose a different path.");
		usage();
		return 1;
	}

	double bounds[6];
	bounds[0] = bounds[2] = bounds[4] = std::numeric_limits<double>::max();
	bounds[1] = bounds[3] = bounds[5] = std::numeric_limits<double>::lowest();

	std::cout << "Computing bounds\n";
	for(const std::string& infile : infiles) {
		geo::pc::PCFile rdr(infile);
		rdr.bounds(bounds);
	}

	double gbounds[6];
	for(int i = 0; i < 6; ++i)
		gbounds[i] = bounds[i];

	int cols, rows;
	std::vector<float> grid;

	std::cout << "Building grid\n";
	if(!buildGrid(infiles, gbounds, resolution, cols, rows, grid))
		return 1;

	geo::pc::PCWriter wtr(outfile, infiles.front());

	std::cout << "Normalizing\n";
	normalize(infiles, wtr, gbounds, resolution, cols, rows, grid);

	return 0;
}
