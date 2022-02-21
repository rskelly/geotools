/*
 * catchment.cpp
 *
 *  Created on: Apr. 24, 2021
 *      Author: rob
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "geo.hpp"
#include "grid.hpp"
#include "util.hpp"

using namespace geo::grid;
using namespace geo::util;

bool flowsinto(int c, int r, unsigned char d, int col, int row) {
	// Return true if the flow from c, r given by d flows into col,row
	if(d == 0) {
		//g_debug("d " << (int) d);
		return true;
	}
	int dc = col - c;
	int dr = row - r;
	unsigned char d1 = 0;
	if(dc == -1 && dr == 1) {
		d1 = 2;
	} else if(dc == -1 && dr == 0) {
		d1 = 1;
	} else if(dc == -1 && dr == -1) {
		d1 = 128;
	} else if(dc == 0 && dr == 1) {
		d1 = 4;
	} else if(dc == 0 && dr == -1) {
		d1 = 64;
	} else if(dc == 1 && dr == 1) {
		d1 = 8;
	} else if(dc == 1 && dr == 0) {
		d1 = 16;
	} else if(dc == 1 && dr == -1) {
		d1 = 32;
	}
	//g_debug("d1 " << d << " " << d1 << " " << (d & d1));
	return (d & d1) != 0;
}

int main(int, char** argv) {

	std::string dirfile(argv[1]);
	std::string outfile(argv[2]);
	std::string seedfile(argv[3]);

	std::vector<std::tuple<int, float, float>> pts;
	{
		std::ifstream seeds(seedfile);
		std::vector<std::string> tmp;
		std::string line;
		while(std::getline(seeds, line)) {
			split(std::back_inserter(tmp), line, ",");
			if(tmp.size() < 3)
				continue;
			pts.emplace_back(atoi(tmp[0].c_str()), atof(tmp[1].c_str()), atof(tmp[2].c_str()));
			tmp.clear();
		}
	}
	if(pts.empty()) {
		std::cerr << "No points given." << std::endl;
		return 1;
	}

	Band<int16_t> dir(dirfile, 0, false, true);
	const GridProps& dprops = dir.props();

	GridProps oprops(dprops);
	oprops.setWritable(true);
	oprops.setDataType(DataType::Int32);
	Band<int> out(outfile, oprops, true);

	int cols = dprops.cols();
	int rows = dprops.rows();
	float nd = dprops.nodata();

	std::vector<int> visited(cols * rows);

	std::list<std::tuple<int, int, int>> q;
	for(const auto& it : pts) {
		g_debug("ID: " << std::get<0>(it));
		q.emplace_back(std::get<0>(it), dprops.toCol(std::get<1>(it)), dprops.toRow(std::get<2>(it)));
		std::fill(visited.begin(), visited.end(), 0);
		while(!q.empty()) {
			const auto& item = q.front();
			q.pop_front();

			int id = std::get<0>(item);
			int col = std::get<1>(item);
			int row = std::get<2>(item);

			out.set(col, row, id);
			
			size_t idx = row * cols + col;
			if(visited[idx] > 1)
				g_debug(id << " " << col << " " << row << " " << visited[idx]);
			if(visited[idx]) 
				continue;
			visited[idx]++;

			for(int r = row - 1; r < row + 2; ++r) {
				for(int c = col - 1; c < col + 2; ++c) {
					if(!(c < 0 || c >= cols - 1 || r < 0 || r >= rows - 1 || (c == col && r == row))) {
						unsigned char d = (unsigned char) dir.get(c, r);
						idx = r * cols + c;
//bool f = flowsinto(c, r, d, col, row);
//						g_debug("flow " << f);
						if(!visited[idx] && flowsinto(c, r, d, col, row))
							q.emplace_back(id, c, r);
					}
				}
			}
		}
	}
}
