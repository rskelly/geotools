/*
 * plotter.hpp
 *
 *  Created on: Jun 10, 2019
 *      Author: rob
 */

#ifndef INCLUDE_PLOTTER_HPP_
#define INCLUDE_PLOTTER_HPP_

#include <list>
#include <tuple>
#include <mutex>

namespace hlrg {
namespace util {

class PlotJob {
public:
	std::string filename, title;
	std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>> items;

	PlotJob();
	PlotJob(const std::string& filename, const std::string& title,
			const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items);

};

class Plotter {
private:
	std::mutex m_pltmtx;		///<! A mutex to protect the plot library.
	std::mutex m_qmtx;			///<! A mutex to protect the plot library.
	std::list<PlotJob> m_queue;	///<! Queue for plot jobs.

public:

	bool hasItems() const;

	void enqueue(const std::string& filename, const std::string& title,
			const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items);

	void plot(const PlotJob& job);

	void process();

};

} // util
} // hlrg

#endif /* INCLUDE_PLOTTER_HPP_ */
