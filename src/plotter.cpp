/*
 * plotter.cpp
 *
 *  Created on: Jun 13, 2019
 *      Author: rob
 */

#include <mutex>
#include <thread>

#include "matplotlibcpp.h"

#include "plotter.hpp"

using namespace geo::plot;

namespace plt = matplotlibcpp;

namespace {

	// Pipe buffer size.
	constexpr size_t BUFFER_SIZE = 1024 * 1024;

	// The singleton instance.
	Plotter* pinst = nullptr;

}

PlotJob::PlotJob() {}

PlotJob::PlotJob(const std::string& filename, const std::string& title,
		const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items) :
	filename(filename), title(title), items(items) {
}


Plotter::Plotter() :
	m_procid(0) {
	if(pinst)
		throw std::runtime_error("There is already an instance of Plotter.");
	pinst = this;
	m_pipefd[0] = 0;
	m_pipefd[1] = 0;
}

Plotter& Plotter::instance() {
	return *pinst;
}

bool Plotter::start() {
	stop();
	if(pipe(m_pipefd))
		return false;
	m_procid = fork();
	if(m_procid == 0) {
		std::cout << "Starting plotter.\n";
		process();
	} else if(m_procid > 0){
		std::cout << "Plotter started.\n";
	} else {
		std::cout << "Failed to start plotter.\n";
		return false;
	}
	return true;
}

bool Plotter::isRunning() {
	return (m_procid > 0 && kill(m_procid, 0) == 0);
}

bool Plotter::stop() {
	if(m_procid > 0) {
		kill(m_procid, SIGKILL);
		close(m_pipefd[0]);
		close(m_pipefd[1]);
		m_procid = 0;
		m_pipefd[0] = 0;
		m_pipefd[1] = 0;
		return true;
	}
	return false;
}

bool Plotter::hasItems() const {
	return !m_queue.empty();
}

void Plotter::enqueue(const std::string& filename, const std::string& title,
		const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items) {
	if(!isRunning()) {
		std::cout << "Not running. Restarting Plotter.\n";
		start();
	}
	std::string output;
	int size;
	{
		std::stringstream ss;
		cereal::BinaryOutputArchive ostr(ss);
		ostr(PlotJob(filename, title, items));
		output = ss.str();
		size = output.size();
	}
	if(write(m_pipefd[1], &size, sizeof(int)) == sizeof(int)) {
		if(write(m_pipefd[1], output.c_str(), output.size()) != size)
			std::cerr << "Failed to enqueue plot job.\n";
	} else {
		std::cerr << "Failed to write size.\n";
		std::cerr << "Failed to enqueue plot job.\n";
	}
}

void Plotter::plot(const PlotJob& job) {
	std::cout << "Plot\n";
	std::lock_guard<std::mutex> lk(m_pltmtx);
	plt::figure(111);
	plt::figure_size(600, 400);
	for(const auto& item : job.items)
		plt::named_plot(std::get<0>(item), std::get<1>(item), std::get<2>(item));
	plt::title(job.title);
	plt::legend();
	plt::save(job.filename);
	plt::close();
}

void parseData(std::vector<char>& buf, int size, std::list<PlotJob>& queue) {
	std::stringstream ss;
	PlotJob job;
	for(int i = 0; i < size; ++i)
		ss << buf[i];
	{
		cereal::BinaryInputArchive istr(ss);
		istr(job);
	}
	queue.push_back(job);
}

void Plotter::process() {
	std::vector<char> buf(BUFFER_SIZE);
	int len, size = 0, offset = 0;
	while(true) {
		if(size == 0) {
			// First, try to read the size.
			len = read(m_pipefd[0], &size, sizeof(int));
			if(len != sizeof(int) || size <= 0)
				size = 0;
		} else {
			// Read the data, up to the amount available.
			len = read(m_pipefd[0], buf.data() + offset, size);
			if(len > 0 && len == size) {
				// If complete, parse the packet.
				parseData(buf, size, m_queue);
				size = 0;
				offset = 0;
			} else if(len > 0 && len < size) {
				// Else set up to read the next part.
				offset += len;
				size -= len;
			}
		}
		while(!m_queue.empty()) {
			PlotJob& job = m_queue.front();
			try {
				plot(job);
			} catch(const std::exception& ex) {
				std::cerr << ex.what() << "\n";
			}
			m_queue.pop_front();
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
	}
}
