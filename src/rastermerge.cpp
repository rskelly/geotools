/*
 * rastermerge.cpp
 *
 *  Created on: Dec 19, 2019
 *      Author: rob
 */

#include <string>
#include <vector>
#include <iostream>
#include <list>
#include <thread>
#include <mutex>

#include <gdal_priv.h>

#include "alglib/interpolation.h"
#include "grid.hpp"
#include "util.hpp"

using namespace geo::grid;
using namespace geo::util;

bool loadRasters(const std::string& outfile, const std::vector<std::string>& files,
		std::vector<int> bands, int count, int resample, Band<float>& data, bool mapped) {

	GridProps props;
	Bounds<double> bounds;
	std::vector<std::unique_ptr<Band<float>>> rasters;
	// Create the list of bands to load.
	for(size_t i = 0; i < (size_t) count; ++i) {
		rasters.emplace_back(new Band<float>(files[i], bands[i] - 1, false, mapped));
		rasters.front()->fixNaNs();
		if(i == 0)
			props = rasters.front()->props();
		bounds.extend(props.bounds());
	}

	// If resampling is set, resize the output raster.
	if(resample > 1) {
		props.setResolution(props.resX() / resample, props.resY() / resample);
		props.setSize(props.cols(), props.rows());
	}
	props.setBounds(bounds);
	props.setWritable(true);
	props.setBands(1);

	// Initialize the output raster.
	data.init(outfile, props, true);
	data.fill(props.nodata());

	// Copy the constituent bands into the output raster.
	for(std::unique_ptr<Band<float>>& band : rasters) {
		const GridProps& p = band->props();
		// If resampling is set, copy pixels one by one.
		for(int row = 0; row < p.rows(); row += resample) {
			for(int col = 0; col < p.cols(); col += resample) {
				float v = band->get(col, row);
				if(v != p.nodata()) { // Leave the output raster's nodata intact if this pixel has nodata.
					int c = props.toCol(p.toX(col));
					int r = props.toRow(p.toY(row));
					data.set(c, r, band->get(col, row));
				}
			}
		}
		/*} else {
			// Copy the whole raster.
			int c = props.toCol(p.toX(0));
			int r = props.toRow(p.toY(0));
			band->writeTo(data, p.cols(), p.rows(), 0, 0, c, r);
		}*/
	}
	return true;
}

// Smooth the raster.
bool smooth(std::vector<bool>& filled, Band<float>& src, Band<float>& dst, int col, int row, int cols, int rows) {
	if(filled[row * cols + col])
		return false;
	float t = 0, w = 0;
	int n = 0;
	double v;
	for(int r = 0; r < rows; ++r) {
		for(int c = 0; c < cols; ++c) {
			if(!std::isnan((v = src.get(c, r)))) {
				float w0 = c == col && r == row ? 1.0 : 1.0 / (std::pow(c - col, 2) + std::pow(r - row, 2));
				t += v * w0;
				w += w0;
				++n;
			}
		}
	}
	if(n && w > 0) {
		dst.set(col, row, t / w);
		filled[row * cols + col] = true;
		return true;
	}
	return false;
}

/**
 * Process inverse distance waiting.
 * \param rowq A queue of row index/block height pairs.
 * \param qmtx A mutex to protect the queue.
 * \param src The source (difference) raster.
 * \param dst The target (output) raster.
 * \param dmtx A mutex to protect the output raster.
 * \param size The search radius in pixels.
 */
void processIDW(std::list<std::pair<int, int>>* rowq, std::mutex* qmtx,
		Band<float>* src, Band<float>* dst, std::mutex* dmtx, int size) {

	// The row index and block height.
	int row, block;

	// The source and destination properties.
	const GridProps& sprops = src->props();
	const GridProps& dprops = dst->props();
	
	// Row data buffer.
	std::vector<float> rowbuf;

	// Cols and rows in the destination raster.
	int dcols = dst->props().cols();
	int drows = dst->props().rows();

	while(!rowq->empty()) {
		{
			// Grab a row/block pair off the queue.
			std::lock_guard<std::mutex> lk(*qmtx);
			if(!rowq->empty()) {
				std::pair<int, int>& item = rowq->front();
				row = item.first;
				block = item.second;
				rowq->pop_front();
				std::cout << "Row " << row << " of " << drows << "\n";
			} else {
				std::this_thread::yield();
				continue;
			}
		}

		// If the block size is too large, shrink it.
		if(row + block >= drows)
			block = drows - row;

		// Resize and fill the row data buffer.
		rowbuf.resize(dcols * block);
		std::fill(rowbuf.begin(), rowbuf.end(), dprops.nodata());

		// Squared radius for validation.
		float maxrad = size * size;
		float nearest;
		// Raster values.
		float v0, v1;
		for(int b = 0; b < block; ++b) {
			for(int col = 0; col < dcols; ++col) {
				// If the destination pixel is invalid skip it.
				if((v1 = dst->get(col, row + b)) == dprops.nodata())
					continue;
				// The sum and weight variables.
				float s = 0;
				float w = 0;
				// If true, exit the loops.
				bool halt = false;
				// Col and row in the difference raster.
				int scol = sprops.toCol(dprops.toX(col));
				int srow = sprops.toRow(dprops.toY(row));
				nearest = geo::maxvalue<float>();
				// Run the circular kernel.
				for(int r = -size; !halt && r < size+ 1; ++r) {
					for(int c = -size; c < size + 1; ++c) {
						// Calculate the kernel pixel.
						int cc = scol + c;
						int rr = srow + r;
						// The distance from the kernel center.
						float d = (float) (c * c + r * r);
						// If the distance is within the max radius, and the pixel is valid, add to sum and weight.
						if(d <= maxrad
								&& sprops.hasCell(cc, rr)
								&& (v0 = src->get(cc, rr)) != sprops.nodata()) {
							if(false) { //d == 0) {
								// If the distance is zero, use the value directly and halt.
								s = v0;
								w = 1;
								halt = true;
								break;
							} else {
								// Add to sum and weight.
								if(d < nearest)
									nearest = d;
								float w0 = 1.0f - (d / maxrad);
								s += v0 * w0;
								w += w0;
							}
						}
					}
				}
				// If the weight is valid, adjust the value.
				if(w > 0) {
					rowbuf[b * dcols + col] = v1 + (s / w) * (1.0f - nearest / maxrad);
				} else {
					rowbuf[b * dcols + col] = v1;
				}
			}
		}

		// Write the block to the output raster.
		std::lock_guard<std::mutex> lk(*dmtx);
		for(int b = 0; b < block; ++b)
			dst->setRow(row + b, rowbuf.data() + (b * dcols));
	}
}

/**
 * Process inverse distance waiting.
 * \param rowq A queue of row index/block height pairs.
 * \param qmtx A mutex to protect the queue.
 * \param src The source (difference) raster.
 * \param dst The target (output) raster.
 * \param dmtx A mutex to protect the output raster.
 * \param size The search radius in pixels.
 */
void processRBF(std::list<std::pair<int, int>>* rowq, std::mutex* qmtx,
		Band<float>* src, Band<float>* dst, std::mutex* dmtx, int size) {

	// The row index and block height.
	int row, block;

	// The source and destination properties.
	const GridProps& sprops = src->props();
	const GridProps& dprops = dst->props();

	// Row data buffer.
	std::vector<float> rowbuf;

	// Cols and rows in the destination raster.
	int dcols = dst->props().cols();
	int drows = dst->props().rows();

	while(!rowq->empty()) {
		{
			// Grab a row/block pair off the queue.
			std::lock_guard<std::mutex> lk(*qmtx);
			if(!rowq->empty()) {
				std::pair<int, int>& item = rowq->front();
				row = item.first;
				block = item.second;
				rowq->pop_front();
				std::cout << "Row " << row << " of " << drows << "\n";
			} else {
				std::this_thread::yield();
				continue;
			}
		}

		// If the block size is too large, shrink it.
		if(row + block >= drows)
			block = drows - row;

		// Resize and fill the row data buffer.
		rowbuf.resize(dcols * block);
		std::fill(rowbuf.begin(), rowbuf.end(), dprops.nodata());

		// Squared radius for validation.
		float maxrad = size * size;
		float nearest;
		// Raster values.
		float v0, v1;
		for(int b = 0; b < block; ++b) {
			for(int col = 0; col < dcols; ++col) {
				// If the destination pixel is invalid skip it.
				if((v1 = dst->get(col, row + b)) == dprops.nodata())
					continue;
				// The sum and weight variables.
				float s = 0;
				float w = 0;
				// If true, exit the loops.
				bool halt = false;
				// Col and row in the difference raster.
				int scol = sprops.toCol(dprops.toX(col));
				int srow = sprops.toRow(dprops.toY(row));
				nearest = geo::maxvalue<float>();
				// Run the circular kernel.
				for(int r = -size; !halt && r < size+ 1; ++r) {
					for(int c = -size; c < size + 1; ++c) {
						// Calculate the kernel pixel.
						int cc = scol + c;
						int rr = srow + r;
						// The distance from the kernel center.
						float d = (float) (c * c + r * r);
						// If the distance is within the max radius, and the pixel is valid, add to sum and weight.
						if(d <= maxrad
								&& sprops.hasCell(cc, rr)
								&& (v0 = src->get(cc, rr)) != sprops.nodata()) {
							if(false) { //d == 0) {
								// If the distance is zero, use the value directly and halt.
								s = v0;
								w = 1;
								halt = true;
								break;
							} else {
								// Add to sum and weight.
								if(d < nearest)
									nearest = d;
								float w0 = 1.0f - (d / maxrad);
								s += v0 * w0;
								w += w0;
							}
						}
					}
				}
				// If the weight is valid, adjust the value.
				if(w > 0) {
					rowbuf[b * dcols + col] = v1 + (s / w) * (1.0f - nearest / maxrad);
				} else {
					rowbuf[b * dcols + col] = v1;
				}
			}
		}

		// Write the block to the output raster.
		std::lock_guard<std::mutex> lk(*dmtx);
		for(int b = 0; b < block; ++b)
			dst->setRow(row + b, rowbuf.data() + (b * dcols));
	}
}

int main(int argc, char** argv) {

	if(argc < 6) {
		std::cerr << "Usage: rastermerge [options] <<anchor file 1> <anchor band 1> [<anchor file 2> <anchor band 2> [...]]> <target file 2> <target band 2> <output file>\n"
				<< " -m <method         avg (shift by average difference), idw or rbf.\n"
				<< " -s <size>          The radius of the window in map units. If a comma-delimited list, will make\n"
				<< "                    multiple passes.\n"
				<< " -r <mult>          A resample multiplier. 1 is no change, 2 doubles the side of a pixel, etc.\n"
				<< "                    If -s is a list, this must be a list of the same length.\n"
				<< " -t <threads>       The number of threads.\n"
				<< " -k <mask> <band>   A mask file. Pixel value 1 is kept.\n"
				<< " -y                 If given, use in-core memory instead of mapped.\n"
				<< " -b <block size>    The height of processed rows. Blocks are processed in parallel. Default 1\n";
		return 1;
	}

	std::vector<std::string> files;
	std::vector<int> bands;
	std::vector<float> radii;
	std::vector<int> resamples;
	int tcount = 4;
	int block = 1;
	std::string outfile;
	std::string maskfile;
	int maskband = 0;
	bool mapped = true;
	std::string method;

	for(int i = 1; i < argc; ++i) {
		std::string arg(argv[i]);
		if(arg == "-s") {
			std::vector<std::string> rad;
			geo::util::split(std::back_inserter(rad), argv[++i], ",");
			for(const std::string& r : rad)
				radii.push_back(atof(r.c_str()));
		} else if(arg == "-m") {
			method = argv[++i];
		} else if(arg == "-r") {
			std::vector<std::string> res;
			geo::util::split(std::back_inserter(res), argv[++i], ",");
			for(const std::string& r : res)
				resamples.push_back(atoi(r.c_str()));
		} else if(arg == "-y") {
			mapped = false;
		} else if(arg == "-b") {
			block = atoi(argv[++i]);
		} else if(arg == "-k") {
			maskfile = argv[++i];
			maskband = atoi(argv[++i]);
		} else if(arg == "-t") {
			tcount = atoi(argv[++i]);
			if(tcount < 1)
				tcount = 1;
		} else {
			if(i < argc - 1) {
				files.push_back(arg);
				bands.push_back(atoi(argv[++i]));
			} else {
				outfile = arg;
			}
		}
	}

	if(radii.empty())
		g_runerr("No radii specified. Need at least one.");
	if(radii.size() != resamples.size())
		g_runerr("The radius and resample lists should be the same length.");

	if(method != "idw")
		radii.resize(1);

	g_debug("Using mapped memory: " << mapped);

	std::string targetFile = files.back();
	int targetBand = bands.back();

	std::string filebase = outfile.substr(0, outfile.find('.'));
	std::cout << filebase << "\n";

	Band<float> output;		// The output raster. Target is loaded, then modified and written.
	GridProps oprops;		// The output raster properties.

	// Load the target/output raster.
	g_debug("Loading target raster.");
	{
		Band<float> target(targetFile, targetBand - 1, false, mapped);
		target.fixNaNs();
		oprops = target.props();
		oprops.setWritable(true);
		oprops.setBands(1);
		output.init(outfile, oprops, mapped);
		target.writeTo(output);
	}

	// Load the anchor rasters(s).
	g_debug("Loading anchor raster.")
	Band<float> anchor;	// The anchor raster.
	loadRasters("/tmp/anchor.tif", files, bands, files.size() - 1, 1, anchor, mapped);
	GridProps aprops = anchor.props();

	// Load the mask raster if available.
	bool hasMask = !maskfile.empty();
	Band<int> mask;		// The mask raster.
	GridProps mprops;	// The mask raster properties.
	if(hasMask) {
		g_debug("Loading mask raster.")
		mask.init(maskfile, maskband - 1, false, mapped);
		mprops = mask.props();
	}

	int valid = 0; //  The number of valid difference pixels.

	for(int i = 0; i < (int)radii.size(); ++i) {
		float radius = radii[i];
		int resample = resamples[i];

		g_debug("Pass " << (i + 1) << ": radius: " << radius << "; resample: " << resample);

		// Initialized the (scaled) difference raster.
		Band<float> rdiff;
		GridProps rprops(oprops);

		if(resample > 1) {
			rprops.setResolution(rprops.resX() * resample, rprops.resY() * resample);
			rprops.setSize(rprops.cols() / resample + 1, rprops.rows() / resample + 1);
		}
		rdiff.init("/tmp/rdiff.tif", rprops, mapped);
		rdiff.fill(rprops.nodata());

		// Calculate the difference between anchor and target rasters.
		g_debug("Calculating difference raster.");
		int rows = oprops.rows();
		int cols = oprops.cols();
		float nd = oprops.nodata();
		float v1, v2, v;
		for(int row1 = 0; row1 < rows; ++row1) {
			if(row1 % 100)
				g_debug("Row " << row1 << " of " << rows);
			for(int col1 = 0; col1 < cols; ++col1) {

				// If the target raster is invalid, skip.
				if((v1 = output.get(col1, row1)) == nd || std::isnan(v1))
					continue;

				// Convert the target cell coords to the anchor and mask cell coords.
				double x = oprops.toX(col1);
				double y = oprops.toY(row1);
				int col2 = aprops.toCol(x);
				int row2 = aprops.toRow(y);
				int mc = mprops.toCol(x);
				int mr = mprops.toRow(y);

				// If the anchor raster has the cell, the mask is valid and the anchor value is valid,
				// don't skip.
				if(!aprops.hasCell(col2, row2)
					|| (v2 = anchor.get(col2, row2)) == aprops.nodata()
					|| std::isnan(v2)
					|| (hasMask && !(mprops.hasCell(mc, mr) && mask.get(mc, mr) == 1)))
					continue;
						
				// Add the difference to the difference raster. If the current value
				// is nodata, set it to zero.
				if((v = rdiff.get(col1 / resample, row1 / resample)) == rprops.nodata())
					v = 0;
				rdiff.set(col1 / resample, row1 / resample, v + v2 - v1);
				++valid;
			}
		}

		// Recalculate the average by the resample factor if it's larger than 1.
		if(resample > 1) {
			g_debug("Normalizing difference raster.")
			rows = rprops.rows();
			cols = rprops.cols();
			valid = 0;
			for(int row = 0; row < rows; ++row) {
				for(int col = 0; col < cols; ++col) {
					if((v = rdiff.get(col, row)) != rprops.nodata()) {
						rdiff.set(col, row, v / (resample * resample));
						++valid;
					}
				}
			}
		}

		{
			// Queue contains start index of row and number of rows.
			g_debug("Building row queue.");
			std::list<std::pair<int, int>> rowq;
			for(int row = 0; row < oprops.rows(); row += block)
				rowq.emplace_back(row, block);

			std::mutex qmtx;
			std::mutex dmtx;
			std::vector<std::thread> threads;

			// Calculate the search radius in columns.
			int size = std::ceil((radius * 2) / std::abs(rdiff.props().resX()));
			// If the resample value is >1, recalculate the cell radius.
			if(resample > 1)
				size = size / resample + 1;
			// The cell radius has to be odd.
			if(size % 2 == 0)
				size++;

			g_debug("Starting processor threads.");
			if(method == "idw") {

				for(int i = 0; i < tcount; ++i)
					threads.emplace_back(processIDW, &rowq, &qmtx, &rdiff, &output, &dmtx, size);

			} else if(method == "rbf") {

				alglib::rbfmodel rm;
				try{
					std::vector<double> pts(valid * 3);
					size_t i = 0;
					double v;
					for(int r = 0; r < rprops.rows(); ++r) {
						for(int c = 0; c < rprops.cols(); ++c) {
							if((v = rdiff.get(c, r)) != rprops.nodata()) {
								pts[i++] = rprops.toX(c);
								pts[i++] = rprops.toY(r);
								pts[i++] = v;
							}
						}
					}

					if(pts.size() > 10000) {
						unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
						std::shuffle(pts.begin(), pts.end(), std::default_random_engine(seed));
						pts.resize(10000);
					}

					alglib::rbfcreate(2, 1, rm);
					alglib::real_2d_array xy;
					xy.setcontent(pts.size() / 3, 3, pts.data());
					alglib::rbfsetpoints(rm, xy, pts.size() / 3);
					alglib::rbfsetalgohierarchical(rm, 10000, 10, 0.001);
					alglib::rbfreport rep;
					alglib::rbfbuildmodel(rm, rep);

				} catch(const std::exception& ex) {
					g_error(ex.what());
					return 1;
				}

				for(int i = 0; i < tcount; ++i)
					threads.emplace_back(processRBF, &rowq, &qmtx, &rdiff, &output, &dmtx, size);

			} else if(method == "avg") {

				tcount = 0;
				const GridStats stats = rdiff.stats();
				for(int r = 0; r < oprops.rows(); ++r) {
					for(int c = 0; c < oprops.cols(); ++c) {
						if((v = output.get(c, r)) != oprops.nodata())
							output.set(c, r, v + stats.mean);
					}
				}

			}

			// Wait for threads to complete.
			for(int i = 0; i < tcount; ++i) {
				if(threads[i].joinable())
					threads[i].join();
			}
		}

	}

	return 0;
}

