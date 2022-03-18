#pragma once

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstring>
#include <cmath>

#include <string>
#include <vector>
#include <array>
#include <list>
#include <tuple>
#include <set>
#include <queue>

#include <limits>
#include <algorithm>
#include <functional>

#include <random>

#include <stdexcept>

#include <armadillo>
#include <igraph/igraph.h>

class epidemic_transport_model
{

public:
	epidemic_transport_model();
	epidemic_transport_model(const std::array<igraph_t,2>&, const std::function<double(const double&)>&, const int&, const int&, const int&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
	epidemic_transport_model(const igraph_t&, const int&, const int&, const int&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
	epidemic_transport_model(const epidemic_transport_model&);

	~epidemic_transport_model();

	epidemic_transport_model& operator=(const epidemic_transport_model&);


	std::list<std::tuple<double, int, int, int, double, double, double>> simulate(const double&, const bool&);
	std::list<std::tuple<double, int, int, int, double, double, double>> simulate(const double&);

	friend std::ostream &operator<<(std::ostream&, const epidemic_transport_model&);

	class event
	{

	public:
		enum class ACTION
		{
			UNSPECIFIED,
			MOVE,
			COMMUNITY_INFECTION,
			TRANSPORT_INFECTION,
			RECOVERY,
			IMMUNITY_LOSS
		};

		event(const double&, const int&, const int&, const event::ACTION&);
		event();
		~event();

		double time;
		int subject;
		int sender;
		event::ACTION action;

	public:
		friend bool operator<(const epidemic_transport_model::event&, const epidemic_transport_model::event&);

		friend std::ostream &operator<<(std::ostream&, const epidemic_transport_model::event&);
	};

private:

	void initialize(const std::array<igraph_t,2>&, const std::function<double(const double&)>&, const int&, const int&, const int&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);

	std::mt19937_64 random_number_engine;

	double simulation_time = std::numeric_limits<double>::infinity();
	std::priority_queue<event> event_queue;

	//igraph_t transport_network;
	std::array<igraph_t,2> transport_networks;
	igraph_t epidemic_network;

	std::function<double(const double&)> transport_network_interpolation_functional;

	// PARAMATERS

	int world_size;
	int community_size;

	int community_network_degree;

	double mobility_rate;
	double community_infection_rate;
	double transport_infection_rate;
	double recovery_rate;
	double immunity_loss_rate;

	int initial_site;
	double initial_prevalence;

	double fractional_exponent;

	// ---------

	std::uniform_real_distribution<double> uniform_distribution = std::uniform_real_distribution<double>(0,1);
	std::exponential_distribution<double> exponential_distribution = std::exponential_distribution<double>(1);

	std::array<std::vector<std::discrete_distribution<int>>,2> transport_transition_distributions;

	std::vector<std::set<int>> site_occupancy;
	std::vector<std::set<int>> community_contacts;

	enum class HEALTH
	{
		S,
		I,
		R
	};

	std::vector<epidemic_transport_model::HEALTH> state_health;
	std::vector<int> state_site;

	std::vector<double> state_site_future_transition_time;
	std::vector<double> state_health_future_recovery_time;
	std::vector<double> state_health_future_immunity_loss_time;

	void infer_transport_transition_distributions(void);

	void initialise_dynamics(void);

	void move_event_handler(const event&);
	void infection_event_handler(const event&);
	void recovery_event_handler(const event&);
	void immunity_loss_event_handler(const event&);

	void infection_spreading(const event&);

	std::list<std::tuple<double, int, int, int, double, double, double>> timeseries = std::list<std::tuple<double, int, int, int, double, double, double>>();

	void update_timeseries(const double&, const bool&);
	void update_timeseries(const double&);

	void print_state(void);
};
