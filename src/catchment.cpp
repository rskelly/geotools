/*
 * catchment.cpp
 *
 *  Created on: Apr. 24, 2021
 *      Author: rob
 */

#include <string>
#include <vector>
#include <iostream>

#include "grid.hpp"
#include "util.hpp"

using namespace geo::grid;
using namespace geo::util;

int main(int argc, char** argv[]) {

	std::string accumfile(argv[1]);
	std::string outfile(argv[2]);
	std::string ptsstr(argv[3]);

	std::vector<std::tuple<int, float, float>> pts;
	{
		std::vector<std::string> tmp;
		geo::util::split(std::back_inserter(tmp), ptsstr, ",");
		for(size_t i = 0; i < tmp.size(); i += 3)
			pts.emplace_back(atoi(tmp[i].c_str()), atof(tmp[i].c_str()), atof(tmp[i + 1].c_str()));
	}
	if(pts.empty()) {
		std::cerr << "No points given." << std::endl;
		return 1;
	}

	Band<float> accum(accumfile, 1, false, true);
	const GridProps& aprops = accum.props();

	GridProps oprops(aprops);
	oprops.setWritable(true);
	oprops.setDataType(DataType::Int32);
	Band<int> out(outfile, oprops, true);

	int cols = aprops.cols();
	int rows = aprops.rows();
	float nd = aprops.nodata();

	std::vector<bool> visited(cols * rows);
	std::fill(visited.begin(), visited.end(), false);

	std::list<std::tuple<int, int, int>> q;
	for(const auto& it : pts) {
		q.emplace_back(std::get<0>(it), aprops.toCol(std::get<1>(it)), aprops.toRow(std::get<2>(it)));

		while(!q.empty()) {
			const auto& item = q.front();
			q.pop_front();

			int id = std::get<0>(item);
			int col = std::get<1>(item);
			int row = std::get<2>(item);
			float a = accum.get(col, row);

			if(a != nd) {
				for(int r = row - 1; r < row + 2; ++r) {
					for(int c = col - 1; c < col + 2; ++c) {
						if(c < 0 || c >= cols - 1 || r < 0 || r >= rows - 1 || c == col || r == row) {
							float a1 = accum.get(c, r);
							size_t idx = r * cols + c;
							if(!visited[idx] && a1 != nd && a1 < a) {
								out.set(c, r, id);
								q.emplace_back(id, c, r);
								visited[idx] = true;
							}
						}
					}
				}
			}
		}
	}


		}
	}


}
