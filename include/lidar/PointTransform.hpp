/*
 * Transform.hpp
 *
 *  Created on: Mar 7, 2019
 *      Author: rob
 */

#ifndef INCLUDE_LIDAR_POINTTRANSFORM_HPP_
#define INCLUDE_LIDAR_POINTTRANSFORM_HPP_

#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

class PointTransform {
private:
	double rotX, rotY, rotZ;
	double transX, transY, transZ;
public:
	PointTransform() :
		PointTransform(0, 0, 0, 0, 0, 0) {
	}

	PointTransform(double rotX, double rotY, double rotZ, double transX, double transY, double transZ) :
		rotX(rotX), rotY(rotY), rotZ(rotZ),
		transX(transX), transY(transY), transZ(transZ) {}

	Vector3d translation() const {
		return Vector3d(transX, transY, transZ);
	}

	Matrix3d rotation() const {
		Matrix3d rot;
		rot = (AngleAxis<double>(rotX, Vector3d::UnitX())
				* AngleAxis<double>(rotY, Vector3d::UnitY())
				* AngleAxis<double>(rotZ, Vector3d::UnitZ())).matrix();
		return rot;
	}

	static std::vector<PointTransform> read(std::istream& input) {
		std::vector<PointTransform> out;
		std::string line;
		while(std::getline(input, line)) {
			std::stringstream ls(line);
			std::vector<std::string> parts(6);
			for(std::string& part : parts)
				std::getline(ls, part, ',');
			out.emplace_back(atof(parts[0].c_str()), atof(parts[1].c_str()), atof(parts[2].c_str()),
					atof(parts[3].c_str()), atof(parts[4].c_str()), atof(parts[5].c_str()));
		}
		return out;
	}

};



#endif /* INCLUDE_LIDAR_POINTTRANSFORM_HPP_ */
