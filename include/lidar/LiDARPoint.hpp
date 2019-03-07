/*
 * LiDARPoint.hpp
 *
 *  Created on: Mar 1, 2019
 *      Author: rob
 */

#ifndef INCLUDE_LIDAR_LIDARPOINT_HPP_
#define INCLUDE_LIDAR_LIDARPOINT_HPP_

#include <Eigen/Core>

/**
 * Abstract class representing a single LiDAR point.
 */
class LiDARPoint {
protected:
	double m_x;
	double m_y;
	double m_z;

	LiDARPoint(double x = 0, double y = 0, double z = 0) :
		m_x(x), m_y(y), m_z(z) {}

	double x() const {
		return m_x;
	}

	void x(double x) {
		m_x = x;
	}

	double y() const {
		return m_y;
	}

	void y(double y) {
		m_y = y;
	}

	double z() const {
		return m_z;
	}

	void z(double z) {
		m_z = z;
	}

	Eigen::Vector3d asVector() const {
		return Eigen::Vector3d(m_x, m_y, m_z);
	}

	void transform(const Eigen::Vector3d& translation, const Eigen::Matrix3d& rotation) {
		Eigen::Vector p = asVector();
		p += translation;
		p *= rotation;
		x(p[0]);
		y(p[1]);
		z(p[2]);
	}

	void transform(const Eigen::Matrix3d& rotation, const Eigen::Vector3d& translation) {
		Eigen::Vector p = asVector();
		p *= rotation;
		p += translation;
		x(p[0]);
		y(p[1]);
		z(p[2]);
	}

	void transform(const Eigen::Vector3d& translation) {
		Eigen::Vector p = asVector();
		p += translation;
		x(p[0]);
		y(p[1]);
		z(p[2]);
	}

	void transform(const Eigen::Matrix3d& rotation) {
		Eigen::Vector p = asVector();
		p *= rotation;
		x(p[0]);
		y(p[1]);
		z(p[2]);
	}

};



#endif /* INCLUDE_LIDAR_LIDARPOINT_HPP_ */
