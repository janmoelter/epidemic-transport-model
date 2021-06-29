#pragma once

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstring>

#include <string>
#include <vector>
#include <array>
#include <list>
#include <tuple>
#include <set>

#include <limits>

#include <algorithm>

#include <queue>

#include <random>

#include <igraph/igraph.h>

enum class epidemic_transport_model_event_action { MOVE, INFECTION, TRANSIENT_INFECTION, RECOVERY, UNSPECIFIED };

class epidemic_transport_model_event
{

public:
	epidemic_transport_model_event(const double&, const int&, const int&, const epidemic_transport_model_event_action&);
	epidemic_transport_model_event();
	~epidemic_transport_model_event();

//private:

	double time;
	int subject;
	int sender;
	epidemic_transport_model_event_action action;

public:
	friend bool operator<(const epidemic_transport_model_event&, const epidemic_transport_model_event&);

	friend std::ostream& operator<<(std::ostream&, const epidemic_transport_model_event&);
};



class epidemic_transport_model
{

public:
	epidemic_transport_model(const std::string&, const int&, const int&, const double&, const double&, const double&, const double&, const double&);
	~epidemic_transport_model();


	std::list<std::tuple<double, int, double>> run(const double&);

private:
	double time_quantum = std::numeric_limits<double>::min();
	double runtime = std::numeric_limits<double>::infinity();

	std::random_device random_number_engine;

	std::priority_queue<epidemic_transport_model_event> event_queue;

	igraph_t transport_network;
	igraph_t epidemic_network;

	// PARAMATER

	double mobility_rate;
	double infection_rate;
	double recovery_rate;
	double coupling;

	// ---------

	double prevalence;

	int community_size;
	int world_size;

	std::exponential_distribution<double> T_transport_transition;
	std::exponential_distribution<double> T_infection;
	std::exponential_distribution<double> T_transient_infection;
	std::exponential_distribution<double> T_recovery;


	std::vector<std::vector<double>> transport_transition_matrix;
	std::vector<std::discrete_distribution<int>> transport_transition_distributions;

	std::vector<std::set<int>> community_adjacents;
	std::vector<std::set<int>> site_occupancy;

	std::vector<double> future__T_move;
	std::vector<double> future__T_recovery;

	void infer_transport_transition_matrix(void);

	void initialise_state(void);


	void move_event_handler(const epidemic_transport_model_event&);
	void infection_event_handler(const epidemic_transport_model_event&);
	void recovery_event_handler(const epidemic_transport_model_event&);

	void infection_spreading(const epidemic_transport_model_event&);

	std::list<std::tuple<double, int, double>> timeseries = std::list<std::tuple<double, int, double>>();

	void update_timeseries(const double&, const bool&);
	void update_timeseries(const double&);

	void print_state(void);

};
