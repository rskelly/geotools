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

namespace geo {
namespace plot {

/**
 * \brief PlotJob contains the information needed to print a plot.
 */
class PlotJob {
public:
	std::string filename, title;															///<! The title of the plot/
	std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>> items;	///<! A list of items which are to be plotted. Each contains a title, a list of abscissae and a list of ordinates.

	/**
	 * \brief Create an empty PlotJob.
	 */
	PlotJob();

	/**
	 * \brief Create a configured PlotJob.
	 *
	 * \param filename The output filename.
	 * \param title The title of the plot.
	 * \param items A list of tuples containing a single plottable item. Each title contains a title, a list of abscissae and a list of ordinates.
	 */
	PlotJob(const std::string& filename, const std::string& title,
			const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items);

	/**
	 * \brief Serialize the PlotJob into or out of a stream.
	 *
	 * \param ar The stream.
	 */
	template <class Archive>
	void serialize(Archive& ar) {
		ar(filename, title, items);
	}

};


/**
 * \brief The Plotter object.
 *
 * The Plotter is responsible for printing PlotJobs. The Plotter runs
 * as a separate process and receives jobs over a pipe.
 */
class Plotter {
private:
	std::mutex m_pltmtx;		///<! A mutex to protect the plot library.
	std::mutex m_qmtx;			///<! A mutex to protect the plot library.
	std::list<PlotJob> m_queue;	///<! Queue for plot jobs.

	int m_pipefd[2];			///<! The pipe file descriptor.
	pid_t m_procid;				///<! The process ID.

	/**
	 * \brief Return true if there are items to plot.
	 *
	 * \return  True if there are items to plot.
	 */
	bool hasItems() const;

	/**
	 * \brief Plot the given PlotJob.
	 *
	 * \param job A PlotJob.
	 */
	void plot(const PlotJob& job);

	/**
	 * \brief Run the processor.
	 */
	void process();

	/**
	 * \brief Return true if the processor is running.
	 *
	 * \return True if the processor is running.
	 */
	bool isRunning();

public:

	/**
	 * \brief Construct the Plotter.
	 */
	Plotter();

	/**
	 * \brief Return the single instance of Plotter.
	 *
	 * \return A reference to the single instance.
	 */
	static Plotter& instance();

	/**
	 * \brief Start the plotter service.
	 *
	 * \return True if the service starts.
	 */
	bool start();

	/**
	 * \brief Stop the plotter service.
	 *
	 * \return True if the service stops.
	 */
	bool stop();

	/**
	 * \brief Add a plot job to the Plotter.
	 *
	 * \param filename The output file.
	 * \param title The plot title.
	 * \param items A vector containing tuples. Each tuple contains a title and one vector each of ordinates and abscissae.
	 */
	void enqueue(const std::string& filename, const std::string& title,
			const std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>>& items);

};

} // util
} // geo

#endif /* INCLUDE_PLOTTER_HPP_ */
