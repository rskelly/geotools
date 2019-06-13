/*
 * plotter.cpp
 *
 *  Created on: Jun 13, 2019
 *      Author: rob
 */

#include <mutex>

#include "matplotlibcpp.h"

#include "plotter.hpp"

using namespace hlrg::util;

namespace {

}

PlotJob::PlotJob() {}
PlotJob::PlotJob(const std::string& filename, const std::string& title,
		const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items) :
	filename(filename), title(title), items(items) {
}


bool Plotter::hasItems() const {
	return !m_queue.empty();
}

void Plotter::enqueue(const std::string& filename, const std::string& title,
		const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items) {
	std::lock_guard<std::mutex> lk(m_qmtx);
	m_queue.emplace_back(filename, title, items);
}

void Plotter::plot(const PlotJob& job) {
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

void Plotter::process() {
	while(!m_queue.empty()) {
		PlotJob& job = m_queue.front();
		plot(job);
		m_queue.pop_front();
	}
}



