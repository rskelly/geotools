/*
 * plotter.hpp
 *
 *  Created on: Jun 10, 2019
 *      Author: rob
 */

#ifndef INCLUDE_PLOTTER_HPP_
#define INCLUDE_PLOTTER_HPP_

#include <mutex>
#include <list>
#include <tuple>

#include "matplotlibcpp.h"

namespace hlrg {

class PlotJob {
public:
	std::string filename, title;
	std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>> items;

	PlotJob() {}
	PlotJob(const std::string& filename, const std::string& title,
			const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items) :
		filename(filename), title(title), items(items) {
	}
};

class Plotter {
private:
	std::mutex m_pltmtx;		///<! A mutex to protect the plot library.
	std::mutex m_qmtx;			///<! A mutex to protect the plot library.
	std::list<PlotJob> m_queue;	///<! Queue for plot jobs.

public:

	bool hasItems() const {
		return !m_queue.empty();
	}

	void queue(const std::string& filename, const std::string& title,
			const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items) {
		std::lock_guard<std::mutex> lk(m_qmtx);
		m_queue.emplace_back(filename, title, items);
	}

	void plot(const PlotJob& job) {
		std::lock_guard<std::mutex> lk(m_pltmtx);
		namespace plt = matplotlibcpp;
		plt::figure(111);
		plt::figure_size(600, 400);
		for(const auto& item : job.items)
			plt::named_plot(std::get<0>(item), std::get<1>(item), std::get<2>(item));
		plt::title(job.title);
		plt::legend();
		plt::save(job.filename);
		plt::close();
	}

	void process() {
		while(!m_queue.empty()) {
			PlotJob& job = m_queue.front();
			plot(job);
			m_queue.pop_front();
		}
	}
};

}

#endif /* INCLUDE_PLOTTER_HPP_ */
