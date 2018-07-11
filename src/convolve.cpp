/*
 * convolve.cpp
 *
 *  Created on: Jul 10, 2018
 *      Author: rob
 */

#include <cmath>
#include <iostream>
#include <vector>

#define PI 3.1415926535

/**
 * Compute the inverse of the Gaussian: retrieve the x that would give
 * the given y. X is returned as an absolute distance from x0.
 *
 * @param sigma The standard deviation of the function.
 * @param y The value of y for which to find x.
 * @param x0 The mean x.
 */
double invGaussian(double sigma, double y, double x0) {
	return sigma * std::sqrt(-2.0 * std::log(y));
}

/**
 * Compute the value of the Gaussian function for a given mean (x0) and
 * x, given sigma and magnitude 1.
 *
 * @param sigma The standard deviation of the function.
 * @param x The current x.
 * @param x0 The mean x.
 */
double gaussian(double sigma, double x, double x0) {
	return std::exp(-0.5 * std::pow((x - x0) / sigma, 2.0));
}


/**
 * The kernel class holds the step values, plus the min
 * and max wavelengths for the range.
 */
class Kernel {
public:
	double min;
	double max;
	std::vector<double> kernel;

	Kernel(double min, double max) :
		min(min), max(max) {}

	Kernel() : Kernel(0, 0) {}

	void resize(int size) {
		kernel.resize(size);
		kernel[0] = 0;
		kernel[size] = 0;
	}

	double& operator[](size_t idx) {
		return kernel[idx];
	}

	/**
	 * Build a 1-dimensional kernel contaning the discretized Gaussian with
	 * size elements. The Gaussian is parameterized by its mean (the bandwidth) and
	 * standard deviation, which is computed from the full width at half magnitude (FWHM)
	 * using the function 2*sqrt(2*ln(2)) * sigma (https://en.wikipedia.org/wiki/Full_width_at_half_maximum).
	 *
	 * @param size The number of elements in the kernel.
	 * @param wl The wavelength in whatever units; must be uniform for all sources.
	 * @param fwhm Full width at half maximum.
	 * @param threshold The minimum value of the function, which will delineate the "ends" at index 0, and size.
	 */
	void build(int size, double wl, double fwhm, double threshold) {
		// Half-width at half maximum.
		double hwhm = fwhm * 0.5;
		// Recover the std. dev from the hwhm.
		double sigma = hwhm / (2.0 * std::sqrt(2.0 * std::log(2.0)));
		// Get the distance from wl for the threshold value.
		double inv = invGaussian(sigma, threshold, wl);
		// Find the bounds for the curve, and the step width.
		min = wl - inv;
		max = wl + inv;
		double step = (max - min) / size;
		// Compute the steps.
		kernel.resize(size + 1);
		for(int i = 1; i < size; ++i)
			kernel[i] = gaussian(sigma, min + step * i, wl);
	}
};


int main(int argc, char** argv) {

	int size = 40;
	Kernel kernel;
	kernel.build(size, 100, 2, 0.001);
	for(double d : kernel.kernel)
		std::cerr << d << "\n";
	return 0;
}
