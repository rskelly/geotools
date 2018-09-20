/*
 * nano_apx_timesync.cpp
 *
 *  Created on: Sep 17, 2018
 *      Author: rob
 */

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "reader.hpp"

using namespace hlrg;

/**
 * Average the values in the given vector.
 */
double _avg(std::vector<uint16_t>& buf) {
	uint32_t sum = 0;
	for(uint16_t v : buf)
		sum += v;
	return (double) sum / buf.size();
}

/**
 * Average the values in the given vector.
 */
double _avg(std::vector<double>& buf) {
	double sum = 0;
	for(double v : buf)
		sum += v;
	return (double) sum / buf.size();
}

std::string _ts2str(long ts) {
	ts /= 1000;
	char buf[32];
	struct tm* dt = localtime(&ts);
	strftime(buf, 32, "%Y/%m/%d %H:%M:%S", dt);
	std::stringstream ss;
	ss << buf << "." << ts % 1000;
	return ss.str();
}

int main() {

	// Proposed algorithm: Since the flame is the largest dataset, we iterate over the rows, using the frame index
	// from the nano to locate the row offset for the corresponding time. Probably should use a b-tree for searching. -- done

	// The frames won't be found in consecutive order because the flame's integration time is lower than the nano's
	// frame period. Get two consecutive frames and attibute the first time to the first half of the intermediate
	// frames and the second time to the others.


	//std::string apx = "/home/rob/Documents/field/2018_sept/hyper_lidar/100128_2018_08_28_19_34_25/imu_gps.txt";
	std::string apx = "/home/rob/Documents/flame/strobe/2/100154_2018_09_19_18_24_49/imu_gps.txt";
	//std::string rast = "/home/rob/Documents/field/2018_sept/hyper_lidar/100128_2018_08_28_19_34_25/raw_7361";
	//std::string rast = "/home/rob/Documents/field/2018_sept/hyper_lidar/100128_2018_08_28_19_34_25/raw_0";
	std::string rast = "/home/rob/Documents/flame/strobe/2/100154_2018_09_19_18_24_49/raw_0";
	///std::string frameIdx = "/home/rob/Documents/field/2018_sept/hyper_lidar/100128_2018_08_28_19_34_25/frameIndex_7361.txt";
	std::string frameIdx = "/home/rob/Documents/flame/strobe/2/100154_2018_09_19_18_24_49/frameIndex_0.txt";
	//std::string frameIdx = "/home/rob/Documents/field/2018_sept/hyper_lidar/100128_2018_08_28_19_34_25/frameIndex_0.txt";
	//std::string flame = "/home/rob/Documents/field/2018_sept/flame/bartier3_AbsoluteIrradiance_13-05-18-419_conv_nano.csv";
	std::string flame = "/home/rob/Documents/flame/strobe/2/strobe2_AbsoluteIrradiance_11-22-08-459_conv_nano.txt";

	int bufCol = 240; // The column in the buffer to use for comparison.
	int cols = 640; // The number of cols in the raster.
	int flameCol = 100; // The flame column to use.

	FrameIndexReader fi(frameIdx);
	IMUGPSReader ir(apx, -7 * 3600 * 1000); // 7 hours DST offset from UTC
	FlameReader fr(flame, 105.7915494 * 1000); // Flame clock offset: +105.7915494s
	Raster raster(rast);

	FlameRow frow0, frow1;
	long gpsTime, actualGpsTime0 = 0, actualGpsTime1 = 0;
	int frame0 = 0, frame1 = 0;
	std::vector<uint16_t> buffer;

	// Get the first frame index.
	int firstIdx;
	fi.getNearestFrame(0, actualGpsTime0, firstIdx);

	std::ofstream out("./tmp2.out", std::ios::out);

	if(fr.next(frow0)) {

		ir.getGPSTime(frow0.utcTime, gpsTime);
		fi.getNearestFrame(gpsTime, actualGpsTime0, frame0);

		while(fr.next(frow1)) {

			ir.getGPSTime(frow1.utcTime, gpsTime);
			fi.getNearestFrame(gpsTime, actualGpsTime1, frame1);

			std::cerr << _ts2str(frow0.utcTime) << ", " << _ts2str(frow1.utcTime) << ", " << _ts2str(actualGpsTime0) << ", " << _ts2str(actualGpsTime1) << "\n";

			if(frame1 > frame0) {
				int half = frame0 + (frame1 - frame0) / 2;
				for(int row = frame0; row < frame1; ++row) {

					// Get the flame row for this raster row.
					FlameRow& frow = row < half ? frow0 : frow1;

					// Read the pixels.
					raster.get(buffer, row - firstIdx);


					// Print a selection of bands from the middle of the row.
					out << actualGpsTime1 << "," << buffer[bufCol] << "," << frow.bands[flameCol] << ", " << _ts2str(actualGpsTime1) << "\n";

					// divide by irrad

					// write to new raster
				}
			}

			frame0 = frame1;
			actualGpsTime0 = actualGpsTime1;
			frow0 = frow1;
		}
	}
}
