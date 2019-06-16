/*
 * plotter.hpp
 *
 *  Created on: Jun 10, 2019
 *      Author: rob
 */

#ifndef INCLUDE_PLOTTER_HPP_
#define INCLUDE_PLOTTER_HPP_

#include <sys/types.h>
#include <unistd.h>

#include <list>
#include <tuple>
#include <mutex>

#include "cereal/cereal.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
#include "cereal/types/tuple.hpp"

namespace hlrg {
namespace util {

class PlotJob {
public:
	std::string filename, title;
	std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>> items;

	PlotJob();
	PlotJob(const std::string& filename, const std::string& title,
			const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items);

	template <class Archive>
	void serialize(Archive& ar) {
		ar(filename, title, items);
	}

};

class Plotter {
private:
	std::mutex m_pltmtx;		///<! A mutex to protect the plot library.
	std::mutex m_qmtx;			///<! A mutex to protect the plot library.
	std::list<PlotJob> m_queue;	///<! Queue for plot jobs.

	int m_pipefd[2];
	pid_t m_procid;

	bool hasItems() const;

	void plot(const PlotJob& job);

	void process();

	bool isRunning();

public:

	Plotter();

	static Plotter& instance();

	/**
	 * Start the plotter service.
	 *
	 * \return True if the service starts.
	 */
	bool start();

	/**
	 * Stop the plotter service.
	 *
	 * \return True if the service stops.
	 */
	bool stop();

	/**
	 * Add a plot job to the Plotter.
	 *
	 * \param filename The output file.
	 * \param title The plot title.
	 * \param items A vector containing tuples. Each tuple contains a title and one vector each of ordinates and abscissae.
	 */
	void enqueue(const std::string& filename, const std::string& title,
			const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items);

};

} // util
} // hlrg

#endif /* INCLUDE_PLOTTER_HPP_ */
