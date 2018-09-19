/*
 * nano_apx_timesync.cpp
 *
 *  Created on: Sep 17, 2018
 *      Author: rob
 */

#include <string>
#include <iostream>
#include <iomanip>

#include "reader.hpp"

double _avg(std::vector<uint16_t>& buf) {
	uint32_t sum = 0;
	for(uint16_t v : buf)
		sum += v;
	return (double) sum / buf.size();
}

double _avg(std::vector<double>& buf) {
	double sum = 0;
	for(double v : buf)
		sum += v;
	return (double) sum / buf.size();
}

int main() {

	// Proposed algorithm: Since the flame is the largest dataset, we iterate over the rows, using the frame index
	// from the nano to locate the row offset for the corresponding time. Probably should use a b-tree for searching. -- done

	// The frames won't be found in consecutive order because the flame's integration time is lower than the nano's
	// frame period. Get two consecutive frames and attibute the first time to the first half of the intermediate
	// frames and the second time to the others.


	std::string apx = "/home/rob/Documents/field/2018_sept/hyper_lidar/100128_2018_08_28_19_34_25/imu_gps.txt";
	std::string rast = "/home/rob/Documents/field/2018_sept/hyper_lidar/100128_2018_08_28_19_34_25/raw_7361";
	//std::string rast = "/home/rob/Documents/field/2018_sept/hyper_lidar/100128_2018_08_28_19_34_25/raw_0";
	std::string frameIdx = "/home/rob/Documents/field/2018_sept/hyper_lidar/100128_2018_08_28_19_34_25/frameIndex_7361.txt";
	//std::string frameIdx = "/home/rob/Documents/field/2018_sept/hyper_lidar/100128_2018_08_28_19_34_25/frameIndex_0.txt";
	std::string flame = "/home/rob/Documents/field/2018_sept/flame/bartier3_AbsoluteIrradiance_13-05-18-419_conv_nano.csv";

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

		ir.getGPSTime(frow0.timestamp, gpsTime);
		fi.getNearestFrame(gpsTime, actualGpsTime0, frame0);

		while(fr.next(frow1)) {

			ir.getGPSTime(frow1.timestamp, gpsTime);
			fi.getNearestFrame(gpsTime, actualGpsTime1, frame1);

			if(frame1 > 0) {
				int half = frame0 + (frame1 - frame0) / 2;
				for(int row = frame0; row < frame1; ++row) {

					// Get the flame row for this raster row.
					FlameRow& frow = row < half ? frow0 : frow1;

					// Read the pixels.
					raster.get(buffer, row - firstIdx);

					// Print a selection of bands from the middle of the row.
					out << actualGpsTime1 << "," << _avg(buffer) << "," << _avg(frow.bands) << "\n";

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
