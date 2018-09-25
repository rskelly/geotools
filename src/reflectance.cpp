/*
 * nano_apx_timesync.cpp
 *
 *  Created on: Sep 17, 2018
 *      Author: rob
 */

#include "../include/reflectance.hpp"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>


#include <QtCore/QObject>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMessageBox>

#include "reader.hpp"
#include "raster.hpp"
#include "ui/reflectance_ui.hpp"


namespace hlrg {

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

void Reflectance::process(const std::string& imuGps, double imuUTCOffset,
		const std::string& rawRad,
		const std::string& frameIdx,
		const std::string& irradConv, double irradUTCOffset,
		const std::string& reflOut,
		bool* running) {

	// Proposed algorithm: Since the flame is the largest dataset, we iterate over the rows, using the frame index
	// from the nano to locate the row offset for the corresponding time. Probably should use a b-tree for searching. -- done

	// The frames won't be found in consecutive order because the flame's integration time is lower than the nano's
	// frame period. Get two consecutive frames and attibute the first time to the first half of the intermediate
	// frames and the second time to the others.

	// Set up data readers and writers.
	FrameIndexReader fi(frameIdx);
	if(!*running) return;
	IMUGPSReader ir(imuGps, imuUTCOffset * 3600000);
	if(!*running) return;
	FlameReader fr(irradConv, irradUTCOffset * 3600000);
	if(!*running) return;
	Raster raster(rawRad);
	if(!*running) return;
	Raster output(reflOut, raster.cols(), raster.rows(), raster.bands(), 0, Float32);

	FlameRow frow0, frow1;
	long gpsTime, actualGpsTime0 = 0, actualGpsTime1 = 0;
	int frame0 = 0, frame1 = 0;
	std::vector<uint16_t> buffer;
	std::vector<double> refl(raster.cols() * raster.bands());

	// Get the first frame index.
	int firstIdx;
	fi.getNearestFrame(0, actualGpsTime0, firstIdx);

	// Set up a writer to log the output for checking.
	//std::ofstream out("./tmp2.out", std::ios::out);

	// Get the first row from the flame data.
	if(*running && fr.next(frow0)) {

		// Get the nearest frame and times for the flame's time.
		ir.getGPSTime(frow0.utcTime, gpsTime);
		fi.getNearestFrame(gpsTime, actualGpsTime0, frame0);

		// Get the next flame row.
		while(*running && fr.next(frow1)) {

			// Get the nearest frame and times for the flame's time.
			ir.getGPSTime(frow1.utcTime, gpsTime);
			fi.getNearestFrame(gpsTime, actualGpsTime1, frame1);

			//std::cerr << _ts2str(frow0.utcTime) << ", " << _ts2str(frow1.utcTime) << ", " << _ts2str(actualGpsTime0) << ", " << _ts2str(actualGpsTime1) << "\n";

			// If the frame has advanced...
			if(*running && frame1 > frame0) {

				// Compute the number of frames that blong to the first time, vs. the second time.
				int half = frame0 + (frame1 - frame0) / 2;

				// Iterate over the rows between the times.
				for(int row = frame0; *running && row < frame1; ++row) {

					// Get the flame row for this raster row.
					FlameRow& frow = row < half ? frow0 : frow1;

					// Read the pixels.
					raster.get(buffer, row - firstIdx);

					// Print a selection of bands from the middle of the row.
					//out << actualGpsTime1 << "," << buffer[bufCol] << "," << frow.bands[flameCol] << ", " << _ts2str(actualGpsTime1) << "\n";

					// For every band and every cell int he buffer, compute the reflectance using the irradiance values.
					for(int b = 0; *running && b < raster.bands(); ++b) {
						for(int c = 0; *running && c < raster.cols(); ++c) {
							size_t i = b * raster.cols() + c;
							refl[i] = buffer[i] / frow.bands[i];
						}
					}

					if(!*running) return;
					// Write to the new raster
					output.write(refl, row);
				}
			}

			// Save the frame/time/etc., before advancing.
			frame0 = frame1;
			actualGpsTime0 = actualGpsTime1;
			frow0 = frow1;
		}
	}
}


int runWithGui(int argc, char **argv) {
	class ReflectanceApp : public QApplication {
	public:
		ReflectanceApp(int& argc, char** argv) : QApplication(argc, argv) {}
		bool notify(QObject *receiver, QEvent *e) {
			try {
				return QApplication::notify(receiver, e);
			} catch(const std::exception &ex) {
				QMessageBox err;
				err.setText("Error");
				err.setInformativeText(QString(ex.what()));
				err.exec();
				return false;
			}
		}
	};

	ReflectanceApp q(argc, argv);
	Reflectance nt;
	ReflectanceForm form(&nt, &q);
	QDialog qform;
	form.setupUi(&qform);
	qform.show();
	return q.exec();
}

} // hlrg


int main(int argc, char** argv) {

	hlrg::runWithGui(argc, argv);

}
