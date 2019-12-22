/*
 * rastermatch.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: rob
 */
#include <grid.hpp>
#include <vector>
#include <sstream>
#include <thread>
#include <iostream>

#include <Eigen/Core>

#include "grid.hpp"

namespace E = Eigen;

using namespace geo::grid;

bool isEdge(std::vector<float>& vec, int dist, float nd) {
	int side = (int) std::sqrt(vec.size());
	int size = (side - 1) / 2;
	int o = 1;
	while(o <= dist) {
		for(int r = size - o; r < size + o + 1; ++r) {
			for(int c = size - o; c < size + o + 1; ++c) {
				if(c < 0 || r < 0 || c >= side || r >= side)
					continue;
				if(vec[r * side + c] == nd)
					return true;
			}
		}
	}
	return false;
}

bool centreValid(std::vector<float>& vec, float nd) {
	int side = (int) std::sqrt(vec.size());
	int size = (side - 1) / 2;
	return vec[size * side + size] != nd;
}

float idw(std::vector<float>& avec, float anodata, std::vector<float>& tvec, float tnodata) {
	int side = (int) std::sqrt(avec.size());
	int size = side / 2;
	float at = 0;
	float aw = 0;
	float tt = 0;
	float tw = 0;
	bool halt = false;
	for(int r = 0; !halt && r < side; ++r) {
		for(int c = 0; c < side; ++c) {
			double av = avec[r * side + c];
			double tv = tvec[r * side + c];
			float d = std::min((float) size, (float) (std::pow((float) c - size, 2.0) + std::pow((float) r - size, 2.0)));
			double w0 = 1.0 - d / (size * size);
			if(av != anodata) {
				at += av * w0;
				aw += w0;
			}
			if(tv != tnodata) {
				tt += tv * w0;
				tw += w0;
			}
		}
	}
	return tw > 0 && aw > 0 ? (at / aw) - (tt / tw) : 0.0;
}

int main(int argc, char** argv) {

	std::vector<std::string> files;
	std::vector<int> bands;
	int size = 251; // edge kernel size

	for(int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-s") {
			size = atoi(argv[++i]);
		} else {
			files.push_back(arg);
			if(i < argc - 1)
				bands.push_back(atoi(argv[++i]));
		}
	}

	Grid<float> anchor(files[0]);
	int aband = bands[0];
	Grid<float> target(files[1]);
	int tband = bands[1];

	const GridProps& aprops = anchor.props();
	const GridProps& tprops = target.props();

	float anodata = aprops.nodata();
	float tnodata = tprops.nodata();

	GridProps oprops = tprops;
	oprops.setBands(1);
	oprops.setWritable(true);
	Grid<float> output(files[2], oprops);
	output.fill(tnodata, 0);

	std::vector<float> avec(size * size);
	std::vector<float> tvec(size * size);

	std::vector<float> gvec(size * size);
	Grid<float>::gaussianWeights(gvec.data(), size, size / 2, 0);

	std::vector<float> outvec(oprops.cols() * oprops.rows());
	std::fill(outvec.begin(), outvec.end(), tnodata);

	for(int row = 0; row < tprops.rows(); ++row) {
		if(row % 10 == 0)
			std::cout << "Row " << row << " of " << tprops.rows() << "\n";
		for(int col = 0; col < tprops.cols(); ++col) {
			double tv = target.get(col, row, tband - 1);
			if(tv != tnodata) {
				target.writeToVector(tvec, col - size / 2, row - size / 2, size, size, tband - 1);
				int ac = aprops.toCol(tprops.toX(col));
				int ar = aprops.toRow(tprops.toY(row));
				double av = anchor.get(ac, ar, aband - 1);
				if(av != anodata) {
					for(int r = 0; r < size; ++r) {
						for(int c = 0; c < size; ++c) {
							double ov = outvec[row * oprops.cols() + col];
							if(ov == tnodata)
								outvec[row * oprops.cols() + col] = 0;
							if(tvec[r * size + c] != tnodata)
								tvec[r * size + c] += (av - tv) * gvec[r * size + c];
						}
					}
					output.readFromVector(tvec, col - size / 2, row - size / 2, size, size, 0);
				}
			}
		}
	}

	return 0;
}


