/*
 * polygonize.cpp
 *
 *  Created on: Apr 13, 2017
 *      Author: rob
 */

#include <grid.hpp>

void usage() {
	std::cerr << "Usage: polygonize <input raster> [output vector]\n"
			<< " -f <format>     The output format. Any of the available\n"
			<< "                 OGR drivers. Default SQLite. If PG, Postgres,\n"
			<< "                 PostgreSQL or PostGIS treats outfile as a connection string.\n"
			<< " -b <band>       The band. Default 1.\n"
			<< " -d              Remove dangles.\n"
			<< " -h              Remove holes.\n"
			<< " -t <t>          The number of threads. Default 3. There are 2 threads for\n"
			<< "                 reading and writing and n threads for processing polygons.\n"
			<< " -l <layer>      The layer name for database formats.\n"
			<< " -i <field>      The field to store the pixel value (id).\n"
			<< " -g <geom>       The name of the geometry field (geom).\n"
			<< " -v <k=v[,k=v]>  A comma-delimited list of fields and values to insert.\n"
			<< " -m <mask image> A mask image to determine the bounds of the input.\n"
			<< " -ids <1,2,3>    A comma-delimited list of IDs to target. If given,\n"
			<< "                 only polygons with these IDs will be extracted.\n";
}

int main(int argc, char** argv) {

	using namespace geo::grid;

	if(argc < 3) {
		usage();
		return 1;
	}

	std::string driver = "SQLite";
	uint16_t band = 1;
	uint16_t threads = 1;
	bool holes = false;
	bool dangles = false;
	std::vector<std::string> args;
	std::string mask;
	std::string layer;
	std::string idField = "id";
	std::string geomField = "geom";
	std::vector<PolygonValue> fieldValues;
	std::vector<int> targetIDs;

	for(int i = 1; i < argc; ++i) {
		std::string v = argv[i];
		if(v == "-f") {
			driver = argv[++i];
		} else if(v == "-l") {
			layer = argv[++i];
		} else if(v == "-i") {
			idField = argv[++i];
		} else if(v == "-g") {
			geomField = argv[++i];
		} else if(v == "-v") {
			std::vector<std::string> pairs;
			std::vector<std::string> kv;
			geo::util::split(std::back_inserter(pairs), argv[++i], ",");
			for(const std::string pair : pairs) {
				geo::util::split(std::back_inserter(kv), pair, "=");
				fieldValues.emplace_back(kv[0], kv[1]);
				kv.clear();
			}
		} else if(v == "ids") {
			std::vector<std::string> ids;
			geo::util::split(std::back_inserter(ids), argv[++i], ",");
			for(const std::string id : ids)
				targetIDs.push_back(atoi(id.c_str()));
		} else if(v == "-b") {
			band = atoi(argv[++i]);
		} else if(v == "-h") {
			holes = true;
		} else if(v == "-d") {
			dangles = true;
		} else if(v == "-t") {
			threads = atoi(argv[++i]);
		} else if(v == "-m") {
			mask = argv[++i];
		} else {
			args.push_back(argv[i]);
		}
	}

	if(args.size() < 2) {
		std::cerr << "Input and output filenames required.\n";
		usage();
		return 1;
	}

	if(args.size() < 3)
		args.push_back("dn");

	if(threads <= 0) {
		std::cerr << "Illegal thread value. Defaulting to 1.\n";
		threads = 1;
	}

	std::vector<std::string> dbdrivers = {"postgres", "postgresql", "pg", "postgis"};
	std::string ldriver = geo::util::lowercase(driver);

	std::string infile = args[0];
	std::string outfile = args[1];

	Band<uint32_t> test(infile, band - 1, false, true);

	bool d3 = false;

	if(std::find(dbdrivers.begin(), dbdrivers.end(), ldriver) == dbdrivers.end()) {
		test.polygonizeToFile(outfile, layer, idField, driver, fieldValues, holes, dangles, d3, targetIDs);
	} else {
		test.polygonizeToTable(outfile, layer, idField, geomField, fieldValues, holes, dangles, d3, targetIDs);
	}
	return 0;
}


