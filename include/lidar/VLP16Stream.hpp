/*
 * VLP16Stream.hpp
 *
 *  Created on: Mar 1, 2019
 *      Author: rob
 */

#ifndef INCLUDE_LIDAR_VLP16STREAM_HPP_
#define INCLUDE_LIDAR_VLP16STREAM_HPP_

#include <fstream>
#include <cmath>
#include <list>
#include <vector>
#include <iostream>

#include "lidar/LiDARStream.hpp"
#include "lidar/LiDARRange.hpp"

#define TO_RAD(x) (x * M_PI / 180.)

/**
 * Lookup table for vertical angles based on the laser position index. There are 16, with indices 0-15.
 */
constexpr float VLP16_VERT_ANGLE[] = {
		TO_RAD(-15.),
		TO_RAD(1.),
		TO_RAD(-13.),
		TO_RAD(-3.),
		TO_RAD(-11.),
		TO_RAD(5.),
		TO_RAD(-9.),
		TO_RAD(7.),
		TO_RAD(-7.),
		TO_RAD(9.),
		TO_RAD(-5.),
		TO_RAD(11.),
		TO_RAD(-3.),
		TO_RAD(13.),
		TO_RAD(-1.),
		TO_RAD(15.)
};

/**
 * Vertical offset of laser from sensor origin in mm.
 */
constexpr float VLP16_VERT_OFFSET[] = {
		11.2,
		-0.7,
		9.7,
		-2.2,
		8.1,
		-3.7,
		6.6,
		-5.1,
		5.1,
		-6.6,
		3.7,
		-8.1,
		2.2,
		-9.7,
		0.7,
		-11.2
};

/**
 * Current VLP-16 can be Strongest or Last in single mode.
 */
enum ReturnMode {
	Strongest = 0x37,
	Last = 0x38,
	Dual = 0x39
};

/**
 * The type of instrument.
 */
enum Instrument {
	VLP16 = 0x22
};


class VLP16Range : public LiDARRange {
protected:
	float m_range;
	float m_time;
	float m_intensity;
	float m_azimuth;
	float m_angle;
	float m_offset;

public:
	VLP16Range() :
		VLP16Range(0, 0, 0, 0, 0, 0) {}

	VLP16Range(float range, float intensity, float time, float azimuth, float angle, float offset) :
		m_range(range), m_time(time), m_intensity(intensity),
		m_azimuth(azimuth), m_angle(angle), m_offset(offset) {}

	void update(float range, float intensity, float time, float azimuth, float angle, float offset) {
		m_range = range;
		m_time = time;
		m_intensity = intensity;
		m_azimuth = azimuth;
		m_angle = angle;
		m_offset = offset;
	}

	void range(float range) {
		m_range = range;
	}

	float range() const {
		return m_range;
	}

	void azimuth(float azimuth) {
		m_azimuth = azimuth;
	}

	float azimuth() const {
		return m_azimuth;
	}

	void time(float time) {
		m_time = time;
	}

	float time() const {
		return m_time;
	}

	bool equals(const VLP16Range& range) const {
		return m_range == range.m_range && m_time == range.m_time && m_intensity == range.m_intensity
				&& m_azimuth == range.m_azimuth && m_angle == range.m_angle && m_offset == range.m_offset;
	}

	void print(std::ostream& out) {
		out << "Range: " << m_range << "; " << m_azimuth << "; " << m_angle << "\n";
	}
};

/**
 *
 * https://velodynelidar.com/docs/notes/63-9276%20Rev%20A%20VLP-16%20Application%20Note%20-%20Packet%20Structure%20&%20Timing%20Definition_Locked.pdf
 */
class VLP16Stream : public LiDARStream {
protected:

	/**
	 * Parse the block flag. Excpetion if it's invalid.
	 */
	bool parseFlag(std::istream& input) {
		unsigned short flag;
		input >> flag;
		return flag == 0xFFEE;
	}

	/**
	 * Parse the 2-byte azimuth value.
	 */
	float parseAzimuth(std::istream& input) {
		unsigned char a, b;
		input >> a;
		input >> b;
		return (float) ((b << 8) | a);
	}

	/**
	 * Return the range in metres.
	 */
	float parseRange(std::istream& input) {
		unsigned char a, b;
		input >> a;
		input >> b;
		return ((b << 8) | a) * 0.002;
	}

	/**
	 * Return the reflectivity.
	 */
	unsigned char parseIntensity(std::istream& input) {
		unsigned char a;
		input >> a;
		return a;
	}

	/**
	 * Parse the packet header.
	 */
	bool parseHeader(std::istream& input) {
		input.seekg(42, std::ios_base::cur); // Skip the header for now.
		return true;
	}

	// ** If only one return in dual mode, adjacent blocks are identical.
	bool parseBlock(std::istream& input, std::vector<VLP16Range>& ranges, int idx, float& azimuth) {

		if(!parseFlag(input))
			return false;

		azimuth = parseAzimuth(input);

		for(int i = 0; i < 16; ++i) {
			float range = parseRange(input);
			unsigned char intensity = parseIntensity(input);
			float angle = VLP16_VERT_ANGLE[i];
			float offset = VLP16_VERT_OFFSET[i];
			ranges[i].update(range, intensity, 0, TO_RAD(azimuth), angle, offset);
		}

		for(int i = 0; i < 16; ++i) {
			float range = parseRange(input);
			unsigned char intensity = parseIntensity(input);
			float angle = VLP16_VERT_ANGLE[i];
			float offset = VLP16_VERT_OFFSET[i];
			ranges[i + 16].update(range, intensity, 0, 0, angle, offset);
		}

		return true;
	}

	ReturnMode parseReturnMode(std::istream& input) {
		unsigned char mode;
		input >> mode;
		return (ReturnMode) mode;
	}

	/**
	 * Return the timestamp in microseconds.
	 */
	float parseTimestamp(std::istream& input) {
		unsigned char a, b, c, d;
		input >> a;
		input >> b;
		input >> c;
		input >> d;
		return (float) ((d << 24) | (c << 16) | (b << 8) | a);
	}

	Instrument parseInstrument(std::istream& input) {
		unsigned char inst;
		input >> inst;
		return (Instrument) inst;
	}

	void fixAzimuths(std::vector<float>& azimuths) {
		for(int i = 0; i < 11; ++i) {
			while(azimuths[i + 1] < azimuths[i])
				azimuths[i + 1] += 360.0;
		}
	}

	bool parsePacket(std::istream& input, std::list<VLP16Range>& final) {

		static std::vector<VLP16Range> ranges(12 * 32);
		static std::vector<VLP16Range> tmp(32);
		static std::vector<float> azimuths(12);
		static float lastAzimuth = -1;
		static int lastCount = 0;

		if(!parseHeader(input))
			return false;

		for(int i = 0; i < 12; ++i) {
			if(!parseBlock(input, ranges, i, azimuths[i]))
				return false;
		}

		float time = parseTimestamp(input);
		ReturnMode rmode = parseReturnMode(input);
		Instrument inst = parseInstrument(input);

		fixAzimuths(azimuths);

		// Interp the azimuths from the previous batch and send them out.
		if(lastCount) {
			float az = azimuths[0];
			while(az < lastAzimuth)
				az += 360.0;
			float interp = lastAzimuth + (az - lastAzimuth) / 2.0;
			for(int i = 0; i < lastCount; ++i) {
				tmp[i].azimuth(interp);
				final.push_back(tmp[i]);
			}
		}

		// Update the times on all ranges.
		for(int i = 0; i < 12 * 32; ++i)
			ranges[i].time(time);

		// Interploate the azimuths on all ranges except the last 16 (32 if dual mode)
		for(int i = 0; i < 11; ++i) {
			float a = azimuths[i] + (azimuths[i + 1] - azimuths[i]) / 2.0;
			for(int j = 16; j < 32; ++j)
				ranges[i * 32 + j].azimuth(a);
		}

		// Send out the finished ranges; keep back the ones that need azimuths interpolated.
		if(rmode == Dual) {
			for(int i = 0; i < 12; i += 2) {
				for(int j = 0; j < 32; ++j) {
					if(i == 10 && j > 15) {
						// Move the last ones to the tmp list to have the azimuth interpolated on next call.
						tmp[j - 16] = ranges[i * 32 + j];
						tmp[j] = ranges[(i + 1) * 32 + j];
					} else {
						const VLP16Range& r0 = ranges[i * 32 + j];
						const VLP16Range& r1 = ranges[(i + 1) * 32 + j];
						final.push_back(r0);
						// In dual mode, if the pairs are equal, they are not separate hits.
						if(!r0.equals(r1))
							final.push_back(r1);
					}
				}
			}
			lastCount = 32;
		} else {
			for(int i = 0; i < 12; ++i) {
				for(int j = 0; j < 32; ++j) {
					if(i == 11 && j > 15) {
						// Move the last ones to the tmp list to have the azimuth interpolated on next call.
						tmp[j - 16] = ranges[i * 32 + j];
					} else {
						final.push_back(ranges[i * 32 + j]);
					}
				}
			}
			lastCount = 16;
		}

		lastAzimuth = azimuths[azimuths.size() - 1];

		return true;
	}

public:

	void load(const std::string& filename, std::list<VLP16Range>& ranges) {
		std::ifstream input(filename);
		while(parsePacket(input, ranges)) {
			for(VLP16Range& r : ranges)
				r.print(std::cout);
		}
	}

};


#endif /* INCLUDE_LIDAR_VLP16STREAM_HPP_ */
