/*
 * pcnorm2.cpp
 *
 *  Created on: May 5, 2019
 *      Author: rob
 */

#include <pdal/StageFactory.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/io/LasReader.hpp>
#include <pdal/io/LasWriter.hpp>
#include <pdal/io/LasHeader.hpp>
#include <pdal/Options.hpp>
#include <pdal/PointRef.hpp>
#include <pdal/io/BufferReader.hpp>
//#include <pdal/StreamPointTable.hpp>

#include <limits>
#include <fstream>

#include "geo.hpp"
#include "util.hpp"
#include "pointcloud.hpp"

using namespace geo::util;

/*
class TestPointTable : public pdal::StreamPointTable {
public:
	TestPointTable(const std::string& file, pdal::PointView& view) :
		StreamPointTable(*view.table().layout(), 0),
		m_count(0),
		m_file(file),
		m_view(view) {
		m_psize = view.table().layout()->pointSize();
		m_buf.resize(m_psize);
	}

protected:
	size_t m_count;
	size_t m_psize;
	std::string m_file;
	std::fstream m_str;
	std::vector<char> m_buf;

	virtual void reset() override {
		m_count = 0;
	}

	virtual char* getPoint(pdal::PointId index) override {
		if(index >= m_count) {
			std::fill(m_buf.begin(), m_buf.end(), 0);
			while(index >= m_count) {
				m_str.seekg(m_count * m_psize, std::ios::seekdir::_S_beg);
				m_str << m_buf;
				++m_count;
			}
			m_str << m_buf;
		}
		m_str.seekg(m_count * m_psize, std::ios::seekdir::_S_beg);
		m_buf << m_str;
		return m_buf.data();
	}

};
*/


/**
 * \brief Build the interpolation grid.
 *
 * Creates a an average grid by adding the point z coordinates to each cell, then
 * dividing by the number of points in the cell.
 *
 * \param infiles A list of LAS files to read.
 * \param bounds The bounds of the grid.
 * \param res The resolution.
 * \param[out] cols The number of columns in the grid.
 * \param[out] rows The number of rows in the grid.
 * \param[out] grid The raster.
 * \param saveGrid Save the raw grid to a local file called pcnorm.grid
 * \return True if nothing screws up.
 */
bool buildGrid(
		const std::vector<std::string>& infiles,
		double* bounds, double res, int& cols, int& rows, std::vector<float>& grid,
		bool saveGrid) {

	// Increase bounds enough to add 2 cells all around.
	bounds[0] -= res * 2;
	bounds[1] += res * 2;
	bounds[2] -= res * 2;
	bounds[3] += res * 2;

	// Get grid size.
	cols = (int) std::ceil((bounds[2] - bounds[0]) / res);
	rows = (int) std::ceil((bounds[3] - bounds[1]) / res);

	std::vector<float> weights(cols * rows);
	size_t count = 0;

	// To compute weighted mean, need weights and values arrays.
	grid.resize(cols * rows);
	std::fill(grid.begin(), grid.end(), 0);
	std::fill(weights.begin(), weights.end(), 0);
	{

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
				row = (int) (py - bounds[1]) / res;

				for(int r = row - 1; r < row + 2; ++r) {
					for(int c = col - 1; c < col + 2; ++c) {
						if(c >= 0 && r >= 0 && c < cols && r < rows) {

							x = bounds[0] + (c * res) + res * 0.5;
							y = bounds[1] + (r * res) + res * 0.5;
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

	// Fill in gaps. Iteratively average the neighbours of invalid
	// cells until there are none left.
	std::cout << "Filling gaps\n";
	int fillCount = 0;
	do {
		fillCount = 0;
		// Loop over the grid.
		for(int r = 0; r < rows; ++r) {
			for(int c = 0; c < cols; ++c) {
				size_t i = r * cols + c;
				if(weights[i] == 0) {
					// If the cell is invalid, loop over the immediate neighbours.
					double w = 0;
					double h = 0;
					for(int rr = r - 1; rr < r + 2; ++rr) {
						for(int cc = c - 1; cc < c + 2; ++c) {
							if(rr >= 0 && rr < rows && cc >= 0 && cc < cols) {
								size_t ii = rr * cols + cc;
								if(weights[ii] > 0) {
									double w0 = std::pow(rr - r, 2) + std::pow(cc - c, 2);
									h += grid[ii] * w0;
									w += w0;
								}
							}
						}
					}
					if(w > 0) {
						grid[i] = h / w;
						++fillCount;
					}
				}
			}
		}
	} while(fillCount);

	if(saveGrid) {
		g_debug("Saving the grid.");
		std::ofstream grd("pcnorm.grid", std::ios::binary);
		for(char g : grid)
			grd << g;
	}
	return true;
}

/**
 * \brief Interpolate the z-coordinate of the given horizontal coordinate.
 *
 * \param x The target x-coordinate.
 * \param y The target y-coordinate.
 * \param x0 The first triangle x-coordinate.
 * \param y0 The first triangle y-coordinate.
 * \param z0 The first triangle z-coordinate.
 * \param x1 The second triangle x-coordinate.
 * \param y1 The second triangle y-coordinate.
 * \param z1 The second triangle z-coordinate.
 * \param x2 The third triangle x-coordinate.
 * \param y2 The third triangle y-coordinate.
 * \param z2 The third triangle z-coordinate.
 */
double bary(double x, double y,
		double x0, double y0, double z0,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	double w0 = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) / ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2));
	double w1 = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) / ((y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2));
	double w2 = 1 - w0 - w1;
	return (w0 * z0) + (w1 * z1) + (w2 * z2);
}

/**
 * \brief Normalize the given LAS files.
 *
 * \param infiles The list of input files.
 * \param wtr The point cloud writer instance.
 * \param bounds The grid bounds.
 * \param res The resolution.
 * \param cols The number of columns in the grid.
 * \param rows The number of rows in the grid.
 * \param grid The grid.
 */
void normalize(const std::vector<std::string>& infiles, const std::string& outfile,
		double* bounds, double res, int cols, int /*rows*/, std::vector<float>& grid) {

	bounds[4] = std::numeric_limits<double>::max();
	bounds[5] = std::numeric_limits<double>::lowest();

	pdal::PointTable table;
	pdal::LasReader reader;
	pdal::PointViewSet viewset;
	pdal::PointViewPtr view;
	pdal::Dimension::IdList dims;
	std::unordered_map<pdal::Dimension::Id, pdal::Dimension::Type> types;
	pdal::LasHeader hdr;

	pdal::Stage* writer;
	pdal::PointTable wtable;

	size_t total_size = 0;
	size_t num = 0;
	for(const std::string& infile : infiles) {
		std::cout << ++num << " of " << infiles.size() << "\n";
		pdal::Options opts;
		opts.add(pdal::Option("filename", infile));
		reader.setOptions(opts);
		reader.prepare(table);
		hdr = reader.header();
		total_size += hdr.pointCount();

		if(num == 1) {
			dims = table.layout()->dims();
			for(pdal::Dimension::Id id : dims) {
				wtable.layout()->registerDim(id, table.layout()->dimType(id));
				const pdal::Dimension::Detail* d = table.layout()->dimDetail(id);
				types[id] = d->type();
			}
		}
	}

	pdal::PointViewPtr wview(new pdal::PointView(wtable));

	num = 0;
	for(const std::string& infile : infiles) {
		std::cout << ++num << " of " << infiles.size() << "\n";

		pdal::Options opts;
		opts.add(pdal::Option("filename", infile));
		reader.setOptions(opts);
		reader.prepare(table);
		hdr = reader.header();
		size_t size = hdr.pointCount();

		viewset = reader.execute(table);
		view = *viewset.begin();
		dims = view->dims();
		view->size();

		std::vector<char> buf(1024);

		for(size_t i = 0; i < size; ++i) {

			double x = view->getFieldAs<double>(pdal::Dimension::Id::X, i);
			double y = view->getFieldAs<double>(pdal::Dimension::Id::Y, i);
			double z = view->getFieldAs<double>(pdal::Dimension::Id::Z, i);

			// Point's home cell.
			int col = (int) (x - bounds[0]) / res;
			int row = (int) (y - bounds[1]) / res;
			// Center of home cell.
			int cx0 = bounds[0] + (col * res) + res * 0.5;
			int cy0 = bounds[1] + (row * res) + res * 0.5;
			int cz0 = grid[row * cols + col];
			// Cell offsets.
			int col0 = x < cx0 ? col - 1 : col + 1;
			int row0 = y < cy0 ? row - 1 : row + 1;
			// Centers of offset cells.
			int cx1 = bounds[0] + (col0 * res) + res * 0.5;
			int cy1 = bounds[1] + (row * res) + res * 0.5;
			int cz1 = grid[row * cols + col0];
			int cx2 = bounds[0] + (col * res) + res * 0.5;
			int cy2 = bounds[1] + (row0 * res) + res * 0.5;
			int cz2 = grid[row0 * cols + col];

			if(cz0 == -9999.0 || cz1 == -9999.0 || cz2 == -9999.0)
				continue;

			// Get the barycentric z
			double nz = bary(x, y, cx0, cy0, cz0, cx1, cy1, cz1, cx2, cy2, cz2);
			// Write the new z value.

			for(pdal::Dimension::Id id : dims) {
				if(id == pdal::Dimension::Id::Z) {
					wview->setField(pdal::Dimension::Id::Z, i, z - nz);
				} else {
					view->getRawField(id, i, (void*) buf.data());
					wview->setField(id, types[id], i, (const void*) buf.data());
				}
			}

			// Adjust bounds.
			if(z < bounds[4]) bounds[4] = z;
			if(z > bounds[5]) bounds[5] = z;
		}
	}

	pdal::BufferReader brdr;
	brdr.addView(wview);

	pdal::StageFactory fact;
	pdal::Options wopts;
	wopts.add("filename", outfile);
	writer = fact.createStage("writers.las");
	writer->setInput(brdr);
	writer->setOptions(wopts);
	writer->prepare(wtable);
	writer->execute(wtable);

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
	bool saveGrid = false;

	int mode = 0;
	for(int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		if(arg == "-f") {
			force = true;
		} else if(arg == "-g") {
			saveGrid = true;
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
		rdr.init();
		rdr.fileBounds(bounds);
	}

	double gbounds[6];
	for(int i = 0; i < 6; ++i)
		gbounds[i] = bounds[i];

	int cols, rows;
	std::vector<float> grid;

	std::cout << "Building grid\n";
	if(!buildGrid(infiles, gbounds, resolution, cols, rows, grid, saveGrid))
		return 1;

	std::cout << "Normalizing\n";
	normalize(infiles, outfile, gbounds, resolution, cols, rows, grid);

	return 0;
}
