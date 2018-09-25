/*
 * nano_timesync.hpp
 *
 *  Created on: Sep 24, 2018
 *      Author: rob
 */

#ifndef INCLUDE_NANO_TIMESYNC_HPP_
#define INCLUDE_NANO_TIMESYNC_HPP_

#include <string>

namespace hlrg {

class NanoTimesync {
public:
	void sync(const std::string& imuGps, double imuUTCOffset,
			const std::string& rawRad,
			const std::string& frameIdx,
			const std::string& irradConv, double irradUTCOffset,
			const std::string& reflOut,
			bool* running);
};


} // hlrg


#endif /* INCLUDE_NANO_TIMESYNC_HPP_ */
