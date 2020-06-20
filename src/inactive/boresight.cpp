/*
 * boresight.cpp
 *
 *  Created on: Mar 1, 2019
 *      Author: rob
 */



#include "lidar/VLP16Stream.hpp"
#include "lidar/PointTransform.hpp"

int main(int argc, char** argv) {

	if(argc < 3) {
		std::cerr << "Usage: boresight <config file> <pcap file [pcap file [pcap file]]> \n";
		return 1;
	}

	std::string config = argv[1];
	std::vector<std::string> pcaps;
	for(int i = 2; i < argc; ++i)
		pcaps.push_back(argv[++i]);

	std::vector<PointTransform> transforms;
	{
		std::ifstream input(config);
		transforms = PointTransform::read(input);
	}

	VLP16Stream vs;
	std::list<VLP16Range> ranges;
	for(const std::string pcap : pcaps) {
		vs.load(pcap, ranges);
		ranges.clear();
	}

}
