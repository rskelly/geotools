/*
 * task.hpp
 *
 *  Created on: Jul 20, 2019
 *      Author: rob
 */

#ifndef INCLUDE_CLOUD_TASK_HPP_
#define INCLUDE_CLOUD_TASK_HPP_

#include <string>

namespace geo {
namespace cloud {

/**
 * \brief A configuration for a Task.
 *
 * Presumably each implementation of Task will have its
 * own implementation of TaskConfig.
 */
class TaskConfig {
public:

	/**
	 * \brief Load the task configuration from a file.
	 *
	 * \param filename The configuration filename.
	 */
	virtual void load(const std::string& filename) = 0;

	/**
	 * \brief Parse the configuration from a string.
	 *
	 * \param config A string containing configuration information.
	 */
	virtual void parse(const std::string& config) = 0;

	virtual ~TaskConfig() {}
};

/**
 * \brief Passed to a Task during execution to collect
 * status updates and other pertinent information
 * that can be conveyed back to the system.
 */
class TaskStatus;

/**
 * \brief Represents a task that can be run as part of a cloud
 * processing job.
 *
 * The Task is able to both estimate the resources it will need to
 * execute, and to perform the execution itself. Tasks may be chained,
 * in which case the resource estimates will be added together.
 *
 * A Task may be the start of a processing chain, or may
 * receive the input from one or more previous tasks. A Task
 * may be the terminus of a chain or continue through one or
 * more subsequent tasks.
 *
 * The outputs of Task(s) are declared in TaskConfig
 * objects and must match the configured inputs to subsequent tasks.
 * For tasks with indeterminate outputs, a modified TaskConfig object
 * can be retrieved through the outputConfig method.
 */
class Task {
public:

	/**
	 * \brief Configure the Task for running.
	 *
	 * The configuration object contains input and output files,
	 * and any parameters required by the task for execution.
	 *
	 * \param config The TaskConfig object.
	 */
	virtual void configure(const TaskConfig& config) = 0;

	/**
	 * \brief Returns the initial configuration for this task.
	 *
	 * \return A reference to the initial TaskConfig for this Task.
	 */
	virtual const TaskConfig& configuration() const = 0;

	/**
	 * \brief Returns the final task configuration for this Task.
	 *
	 * If the outputs of the Task are indeterminate, the configuration
	 * may have been modified to reflect the outputs. If not,
	 * this method returns the same object as configuration().
	 *
	 * \return A TaskConfig object containing the final configuration including indeterminate outputs.
	 */
	virtual const TaskConfig& outputConfiguration() const = 0;

	/**
	 * \brief Estimate the amount of memory required by the task for execution.
	 *
	 * This is to be considered RAM, and will assume that an indeterminate
	 * amount of RAM is available for execution.
	 *
	 * \return The amount of memory requried for the task to complete.
	 */
	virtual size_t estimateWorkingMemory() const = 0;

	/**
	 * \brief Estimate the amount of 'convertible' memory required by the Task for execution, in addition to working memory.
	 *
	 * Convertible memory is space that would be reserved in RAM, but
	 * can be transferred to disk if necessary. There are two scenarios
	 * where this could be useful: a) when a process has indeterminate needs
	 * and may need to "flip" over to disk-backed storage; b) when the
	 * working memory requirement simply exceeds the maximum available
	 * physical storage and must be forced onto disk.
	 *
	 * If all of the working memory of the Task is convertible, the estimate of
	 * working memory should return zero. In general, this estimate is in addition
	 * to physical working memory.
	 *
	 * \return The amount of convertible working memory.
	 */
	virtual size_t esitmateConvertibleWorkingMemory() const = 0;
	/**
	 * \brief Estimate the on-disk size of the outputs from the task.
	 *
	 * \return An estimate the on-disk size of the outputs from the task.
	 */
	virtual size_t esitmateOutputSize() const = 0;

	/**
	 * \brief An estimate of the time cost of this task.
	 *
	 * The time cost isn't a wall-clock time, but a score which gives
	 * some indication of the time cost of one task relative to another.
	 *
	 * The actual time that a task will take will have to be calculated from the
	 * score based on empirical observation of previous runs of the task.
	 *
	 * \return An estimate of the time cost of the task.
	 */
	virtual float estimateTime() const = 0;

	/**
	 * \brief Execute the task.
	 *
	 * \param status A TaskStatus object to keep track of the state of the task.
	 */
	virtual void execute(TaskStatus& status) = 0;

	virtual ~Task() {}

};


/**
 * \brief Used to track the status of a Job.
 */
class JobStatus;


/**
 * \brief A Job is a single processing job that may contain one or more Tasks.
 *
 * A Job may contain a single Task, a linear chain, or a graph. In the latter
 * case, multiple tasks may produce inputs to a single task, and
 * a single task may provide inputs to multiple tasks.
 *
 * When multiple Tasks are added as predecessors to a Task, their
 * inputs may not all be used by the successor, but all the inputs
 * required by the successor must be produced by the time the
 * task executes, unless they exist at start.
 *
 * Similarly, a Task with multiple successors must provide all inputs
 * to the successors that do not already exist.
 *
 * Tasks do not provide inputs directly to successors. They produce output
 * according to their configuration; successors take those outputs
 * as inputs according to their configuration.
 */
class Job {
public:

	/**
	 * \brief Add a Task before the given task.
	 *
	 * \param task A Task.
	 * \param before The Task which will follow the added Task.
	 */
	void addTaskBefore(const Task& task, const Task& before);

	/**
	 * \brief Add a Task after the given task.
	 *
	 * \param task A Task.
	 * \param after The Task which will preceed the added Task.
	 */
	void addTaskAfter(const Task& task, const Task& after);

	/**
	 * \brief Add a Task.
	 *
	 * This Task will be executed as a starting point, and may have
	 * successors.
	 *
	 * \param task A Task.
	 */
	void addTask(const Task& task);

	/**
	 * \brief Estimate the amount of memory required by the Job for execution.
	 *
	 * This is to be considered RAM, and will assume that an indeterminate
	 * amount of RAM is available for execution.
	 *
	 * \return The amount of memory requried for the Job to complete.
	 */
	size_t estimateWorkingMemory() const;

	/**
	 * \brief Estimate the amount of 'convertible' memory required by the Job for execution, in addition to working memory.
	 *
	 * Convertible memory is space that would be reserved in RAM, but
	 * can be transferred to disk if necessary. There are two scenarios
	 * where this could be useful: a) when a process has indeterminate needs
	 * and may need to "flip" over to disk-backed storage; b) when the
	 * working memory requirement simply exceeds the maximum available
	 * physical storage and must be forced onto disk.
	 *
	 * If all of the working memory of the Job is convertible, the estimate of
	 * working memory should return zero. In general, this estimate is in addition
	 * to physical working memory.
	 *
	 * \return The amount of convertible working memory.
	 */
	size_t esitmateConvertibleWorkingMemory() const;

	/**
	 * \brief Estimate the on-disk size of the outputs from the Job.
	 *
	 * \return An estimate the on-disk size of the outputs from the Job.
	 */
	size_t esitmateOutputSize() const;

	/**
	 * \brief An estimate of the time cost of this Job.
	 *
	 * The time cost isn't a wall-clock time, but a score which gives
	 * some indication of the time cost of one Job relative to another.
	 *
	 * The actual time that a Job will take will have to be calculated from the
	 * score based on empirical observation of previous runs of the Tasks that
	 * comprise the Job.
	 *
	 * \return An estimate of the time cost of the Job.
	 */
	float estimateTime() const;

	/**
	 * \brief Execute the Job.
	 *
	 * \param status A JobStatus object to keep track of the state of the Job.
	 */
	void execute(JobStatus& status);
};

} // cloud
} // geo




#endif /* INCLUDE_CLOUD_TASK_HPP_ */
