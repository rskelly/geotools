/*
 * nano_timesync.hpp
 *
 *  Created on: Sep 24, 2018
 *      Author: rob
 */

#ifndef INCLUDE_REFLECTANCE_HPP_
#define INCLUDE_REFLECTANCE_HPP_

#include <string>

namespace hlrg {

class Reflectance;

class ReflectanceListener {
public:
	virtual void started(Reflectance*) = 0;
	virtual void update(Reflectance*) = 0;
	virtual void stopped(Reflectance*) = 0;
	virtual void finished(Reflectance*) = 0;
	virtual ~ReflectanceListener() {}
};

class Reflectance {
private:
	int m_step;			///<! The current processing step.
	int m_numSteps;		///<! The number of processing steps to complete.

public:
	/**
	 * Processes the radiance image and convolved irradiance spectra to produce a
	 * reflectance image with the same characteristics as the radiance image.
	 *
	 * @param imuGps The imu_gps.txt file produced by the APX-15. This provides a mapping between the GPS
	 * timestamp and the UTC date/time.
	 * @param imuUTCOffset This provides a way of offsetting the UTC timestamp from the APX. Usually not needed. Decimal hours.
	 * @param rawRad The raw radiance image.
	 * @param frameIdx The frame index file produced by the Nano Hyperspec.
	 * @param irradConv The irradiance spectra convolved to correspond to the band map of the Nano.
	 * @param irradUTCOffset A time offset to convert the timestamps in the irradiance spectrometer to match those from the APX-15.
	 * @param reflOut An output image for the reflectance.
	 * @param running Method will continue so long as this is set to true or until completion.
	 */
	void run(ReflectanceListener& listener,
			const std::string& imuGps, double imuUTCOffset,
			const std::string& rawRad,
			const std::string& frameIdx,
			const std::string& irradConv, double irradUTCOffset,
			const std::string& reflOut,
			bool& running);

	double progress() const;

};


} // hlrg


#endif /* INCLUDE_REFLECTANCE_HPP_ */
