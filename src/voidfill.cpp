/*
 * voidfill.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include <thread>
#include <vector>
#include <string>

#include "grid.hpp"

using namespace geo::grid;

void usage() {
	std::cout << "Usage: voidfill [options] <input raster> <output raster>\n"
			<< " -b  <band>       The band. Default 1.\n"
			<< " -m  <area>       Maximum area to fill. Square map units.\n"
			<< " -e               Fill pixels connected to edges. Default off.\n"
			<< " -d  <mode>       Mode: 0=min, 1=mean, 2=median, 3=max. Default 0.\n"
			<< " -n  <n>          The radius of the alpha disc. Implies use of concave hull. Overrides -e.\n"
			<< " -f               Force the overwriting of existing files.\n";
}

/**
 * \brief Fills the voids in a raster using the given mask.
 *
 * \param mask The mask. Pixels with value 1 are filled; 0 not.
 * \param inrast The input raster.
 * \param outrast The output raster.
 * \param col0 The first column of the subset of the raster to be filled.
 * \param row0 The first row of the subset of the raster to be filled.
 * \param col1 The last column (inclusive) of the subset of the raster to be filled.
 * \param row1 The last row (inclusive) of the subset of the raster to be filled.
 */
void fillVoids(Band<uint8_t>& mask, Band<float>& inrast, Band<float>& outrast,
		int col0, int row0, int col1, int row1, int mode) {

	const GridProps& props = inrast.props();
	float nd = props.nodata();
	float v, s;
	int ct = 0;
	std::vector<float> v0;

	// Initialize the values.
	switch(mode) {
	case 0: s = geo::maxvalue<float>(); break;
	case 3: s = geo::minvalue<float>(); break;
	default: s = 0; break;
	}

	// Find the edge pixel values.
	for(int r = row0; r <= row1; ++r) {
		for(int c = col0; c <= col1; ++c) {
			if(mask.get(c, r) != 2)
				continue;
			// Check the pixels around the null pixel for values.
			for(int rr = r - 1; rr < r + 2; ++rr) {
				for(int cc = c - 1; cc < c + 2; ++cc) {
					if(props.hasCell(cc, rr) && (v = inrast.get(cc, rr)) != nd) {
						switch(mode) {
						case 0:
							if(v < s) s = v;
							break;
						case 3:
							if(v > s) s = v;
							break;
						case 2:
							v0.push_back(v);
							break;
						default:
							s += v;
							break;
						}
						++ct;
					}
				}
			}
		}
	}

	float m = nd;

	// Calculate the output value.
	// TODO: Add a spline, IDW (etc.) interp method.
	if(ct) {
		switch(mode) {
		case 0:
		case 3:
			m = s;
			break;
		case 2:
			if(v0.size() % 2 == 0) {
				std::sort(v0.begin(), v0.end());
				m = (v0[v0.size() / 2 - 1] + v0[v0.size() / 2]) / 2.0;
			} else {
				m = v0[v0.size() / 2];
			}
			break;
		default:
			m = s / ct;
			break;
		}
	}

	// Write the new value to the null region, set mask to 3.
	for(int r = row0; r <= row1; ++r) {
		for(int c = col0; c <= col1; ++c) {
			if(mask.get(c, r) == 2) {
				outrast.set(c, r, m);
				mask.set(c, r, 3);
			}
		}
	}
}

/**
 * \brief Set all edge-contacting null pixel regions to 2 in the mask.
 *
 * Iterates over the edge pixels of the raster. Any pixels that are null
 * are flood-filled to transfer the shape of the null region to the mask.
 *
 * \param rast The raster containing voids.
 * \param mask The mask.
 */
void edgeMask(Band<float>& rast, Band<uint8_t>& mask) {
	float nd = rast.props().nodata();
	int cols = rast.props().cols();
	int rows = rast.props().rows();
	TargetFillOperator<float, uint8_t> op1(&rast, 1, &mask, 1, nd, 2); // Mark for filling (nd --> 2)
	int minc, minr, maxc, maxr, area;
	mask.fill(0);
	for(int c = 0; c < cols; ++c) {
		if(mask.get(c, 0) != 2)
			rast.floodFill(c, 0, op1, false, &minc, &minr, &maxc, &maxr, &area);
		if(mask.get(c, rows - 1) != 2)
			rast.floodFill(c, rows - 1, op1, false, &minc, &minr, &maxc, &maxr, &area);
	}
	for(int r = 0; r < rows; ++r) {
		if(mask.get(0, r) != 2)
			rast.floodFill(0, r, op1, false, &minc, &minr, &maxc, &maxr, &area);
		if(mask.get(cols - 1, r) != 2)
			rast.floodFill(cols - 1, r, op1, false, &minc, &minr, &maxc, &maxr, &area);
	}
}

/**
 * \brief Set all non-edge, null pixel regions to 1 in the mask; zero everywhere else.
 *
 * \param rast The raster conatining voids.
 * \param mask The mask.
 */
void ndMask(Band<float>& rast, Band<uint8_t>& mask) {

	edgeMask(rast, mask);

	float nd = rast.props().nodata();
	int cols = rast.props().cols();
	int rows = rast.props().rows();
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {
			if(mask.get(c, r) != 2 && rast.get(c, r) == nd) {
				mask.set(c, r, 1);
			} else {
				mask.set(c, r, 0);
			}
		}
	}
}

/**
 * \brief Process the geometric mask in parallel.
 *
 * \param mask The mask containing pixels to process.
 * \param filled The grid of pixels to be filled.
 * \param kernel The round kernel.
 * \param mtx A mutex to protect the row and filled pointers.
 * \param cols The number of columns.
 * \param rows The number of rows.
 */
void process(const std::vector<bool>* mask, std::vector<bool>* filled,
		const std::vector<std::pair<int, int>>* kernel,
		std::mutex* mtx, int* row, int cols, int rows) {

	// To track which pixels are filled.
	std::vector<bool> filled_(cols * rows);
	std::fill(filled_.begin(), filled_.end(), false);

	// To track which rows should be updated.
	std::vector<int> rowsfilled;

	// Fill the mask with 3 outside the alpha shape.
	int statusStep = std::max(1, rows / 10);
	int r, cc, rr;

	while(true) {
		{
			std::lock_guard<std::mutex> lk(*mtx);
			r = (*row)++;
			if(r >= rows)
				break;
		}
		rowsfilled.push_back(r);
		if(r % statusStep == 0)
			g_debug("Building hull mask. Row: " << r << " of " << rows);
		for(int c = 0; c < cols; ++c) {
			if(mask->at(r * cols + c)) {
				// If all pixels in the fast mask under the disk are true,
				// fill the disk in the mask raster.
				bool fill = true;
				for(const auto& it : *kernel) {
					cc = c + it.first;
					rr = r + it.second;
					if(!(cc < 0 || rr < 0 || cc >= cols || rr >= rows) && !mask->at(rr * cols + cc)) {
						fill = false;
						break;
					}
				}
				if(fill) {
					for(const auto& it : *kernel) {
						cc = c + it.first;
						rr = r + it.second;
						if(!(cc < 0 || rr < 0 || cc >= cols || rr >= rows) && !filled_.at(rr * cols + cc))
							filled_.at(rr * cols + cc) = true;
					}
				}
			}
		}
	}
	std::lock_guard<std::mutex> lk(*mtx);
	for(const int& r : rowsfilled) {
		for(int c = 0; c < cols; ++c) {
			if(filled_.at(r * cols + c))
				filled->at(r * cols + c) = true;
		}
	}
}

/**
 * \brief Set all non-edge, null pixel regions to 1 in the mask.
 *
 * Boundaries are set according to the "alpha" parameter; zero everywhere else.
 *
 * \param rast The raster containing voids.
 * \param mask The mask.
 * \param alpha The "alpha" parameter (just the radius of the "disk").
 */
void hullMask(Band<float>& rast, Band<uint8_t>& mask, double alpha) {

	// Mask the raster's void edges.
	edgeMask(rast, mask);

	float nd = rast.props().nodata();
	int cols = rast.props().cols();
	int rows = rast.props().rows();
	int rad = (int) geo::sq(std::ceil(std::abs(alpha / rast.props().resX())));

	// A fast mask for edge-connected null pixels.
	std::vector<bool> msk(cols * rows);
	std::fill(msk.begin(), msk.end(), false);
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c)
			msk[r * cols + c] = mask.get(c, r) == 2 && rast.get(c, r) == nd;
	}

	// Keeps track of which mask pixels should be filled. (Updated by threads.)
	std::vector<bool> filled(cols * rows);
	std::fill(filled.begin(), filled.end(), false);

	// Build the round kernel.
	std::vector<std::pair<int, int>> k;
	for(int rr = -rad; rr < rad + 1; ++rr) {
		for(int cc = -rad; cc < rad + 1; ++cc) {
			if(geo::sq(rr) + geo::sq(cc) <= rad)
				k.emplace_back(cc, rr);
		}
	}

	// Run the hull generation in parallel.
	int t = std::thread::hardware_concurrency();
	int row = 0;
	std::mutex mtx;
	std::vector<std::thread> threads;
	for(int i = 0; i < t; ++i)
		threads.emplace_back(process, &msk, &filled, &k, &mtx, &row, cols, rows);
	for(int i = 0; i < t; ++i) {
		if(threads[i].joinable())
			threads[i].join();
	}

	// Set hull and nd pixels to 1, valid and edge pixels to 0.
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c)
			mask.set(c, r, filled[r * cols + c] || rast.get(c, r) != nd ? 0 : 1);
	}
	mask.flush();
}

int main(int argc, char** argv) {

	if(argc < 3) {
		usage();
		return 1;
	}

	std::string infile;
	std::string outfile;
	std::string maskfile;
	bool saveMask = false;
	int band = 1;
	float maxarea = geo::maxvalue<float>();
	int mode = 0;
	bool noEdges = true;
	int state = 0;
	float n = 0;
	bool force = false;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-b") {
			band = atoi(argv[++i]);
		} else if(v == "-m") {
			maxarea = atof(argv[++i]);
		} else if(v == "-d") {
			mode = atoi(argv[++i]);
		} else if(v == "-n") {
			n = std::abs(atof(argv[++i]));
		} else if(v == "-f") {
			force = true;
		} else if(v == "-e") {
			noEdges = false;
		} else if(state == 0) {
			infile = v;
			++state;
		} else if(state == 1) {
			outfile = v;
			++state;
		} else if(state == 2) {
			maskfile = v;
			saveMask = true;
			++state;
		}
	}

	if(band - 1 < 0) {
		g_error("Illegal band number: " << band);
		usage();
		return 1;
	}

	if(infile.empty() || outfile.empty()) {
		g_error("Input and output filenames required.");
		usage();
		return 1;
	}

	if(!isfile(infile)) {
		g_error("The input file " << infile << " isn't a valid input file.");
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

	// Load the input raster.
	Band<float> inrast(infile, band - 1 , false, true);

	// Configure the output raster and write the input to it.
	Band<float> outrast(outfile, inrast.props(), true);
	inrast.writeTo(outrast);

	// Configure the mask.
	GridProps mprops(inrast.props());
	mprops.setWritable(true);
	mprops.setBands(1);
	mprops.setDataType(DataType::Byte);
	Band<uint8_t> mask;
	if(saveMask) {
		mask.init(maskfile, mprops, true);
	} else {
		mask.init(mprops, true);
	}
	mask.fill(0);

	int cols = mprops.cols();
	int rows = mprops.rows();

	// Configure operators for flood fill. The first fills 1 pixels with 2,
	// the second fills 2 pixels with 3.
	geo::grid::TargetFillOperator<uint8_t, uint8_t> op1(&mask, 1, &mask, 1, 1, 2); // Mark for filling (1 --> 2)
	geo::grid::TargetFillOperator<uint8_t, uint8_t> op2(&mask, 1, &mask, 1, 2, 3); // Mark for filling (2 --> 3)
	int cmin = 0, cmax = 0, rmin = 0, rmax = 0, area = 0;

	if(n > 0){
		// Build concave hull to produce mask.
		g_debug("Building concave hull mask.");
		hullMask(inrast, mask, n);
	} else if(noEdges){
		// Mask the edge-contacting regions.
		g_debug("Building edge mask.");
		ndMask(inrast, mask);
	} else {
		g_debug("Filling to edges.");
	}

	// Calculate the maximum void fill area (pixels)
	// from the resolution (map units).
	int maxpxarea;
	double q = maxarea / geo::sq(mprops.resX());
	if(q > (double) geo::maxvalue<int>()) {
		maxpxarea = geo::maxvalue<int>();
	} else {
		maxpxarea = (int) q;
	}
	if(maxpxarea < 1) maxpxarea = 1;

	// Do it!
	int statusStep = std::max(1, rows / 10);
	for(int row = 0; row < rows; ++row) {
		if(row % statusStep == 0)
			g_debug("Filling voids. Row: " << row << " of " << rows);
		for(int col = 0; col < cols; ++col) {

			// Only hit pixels marked with 1.
			if(mask.get(col, row) != 1)
				continue;

			// Fill the target region with 2.
			Band<uint8_t>::floodFill(col, row, op1, false, &cmin, &rmin, &cmax, &rmax, &area);

			if(area >= maxpxarea) {
				// The filled area was too large. Ignore it by setting it to 3.
				Band<uint8_t>::floodFill(col, row, op2, false, &cmin, &rmin, &cmax, &rmax, &area);
			} else {
				// Fill the voids.
				fillVoids(mask, inrast, outrast, cmin, rmin, cmax, rmax, mode);
			}
		}
	}

	return 0;
}


