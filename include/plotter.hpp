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

#include "matplotlibcpp.h"

namespace hlrg {

class PlotJob {
public:
	std::string filename, title;
	std::vector<double> x, y, regx, regy;
	PlotJob() {}
	PlotJob(const std::string& filename, const std::string& title,
			const std::vector<double>& x, const std::vector<double>& y,
			const std::vector<double>& regx, const std::vector<double>& regy) :
		filename(filename), title(title), x(x), y(y), regx(x), regy(y) {
	}
};

class Plotter {
private:
	std::mutex m_pltmtx;		///<! A mutex to protect the plot library.
	std::mutex m_qmtx;		///<! A mutex to protect the plot library.
	std::list<PlotJob> m_queue;

public:

	bool hasItems() const {
		return !m_queue.empty();
	}

	void queue(const std::string& filename, const std::string& title,
			const std::vector<double>& x, const std::vector<double>& y,
			const std::vector<double>& regx, const std::vector<double>& regy) {
		std::lock_guard<std::mutex> lk(m_qmtx);
		m_queue.emplace_back(filename, title, x, y, regx, regy);
	}

	void plot(const PlotJob& job) {
		std::lock_guard<std::mutex> lk(m_pltmtx);
		namespace plt = matplotlibcpp;
		plt::figure(111);
		plt::figure_size(600, 400);
		plt::named_plot("Normalized, Continuum Removed", job.x, job.y);
		plt::named_plot("Regression", job.regx, job.regy);
		plt::title(job.title);
		plt::legend();
		plt::save(job.filename);
		plt::close();
	}

	void step() {
		PlotJob job;
		bool hasJob = false;
		while(!m_queue.empty()) {
			{
				std::lock_guard<std::mutex> lk(m_qmtx);
				if(!m_queue.empty()) {
					job = std::move(m_queue.front());
					hasJob = true;
					m_queue.pop_front();
				}
			}
			if(hasJob)
				plot(job);
			hasJob = false;
		}
	}
};

}

#endif /* INCLUDE_PLOTTER_HPP_ */
