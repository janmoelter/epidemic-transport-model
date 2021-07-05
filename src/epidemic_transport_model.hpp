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

#include <random>

#include <igraph/igraph.h>

class epidemic_transport_model
{

public:
	epidemic_transport_model(const igraph_t *, const int &, const int &, const double &, const double &, const double &, const double &, const double &, const double &);
	~epidemic_transport_model();


	std::list<std::tuple<double, int, int, int, double, double, double>> simulate(const double &, const bool &);
	std::list<std::tuple<double, int, int, int, double, double, double>> simulate(const double &);

	friend std::ostream &operator<<(std::ostream &, const epidemic_transport_model &);

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

		event(const double &, const int &, const int &, const event::ACTION &);
		event();
		~event();

		double time;
		int subject;
		int sender;
		event::ACTION action;

	public:
		friend bool operator<(const epidemic_transport_model::event &, const epidemic_transport_model::event &);

		friend std::ostream &operator<<(std::ostream &, const epidemic_transport_model::event &);
	};

private:
	double simulation_time = std::numeric_limits<double>::infinity();

	std::random_device random_number_engine;

	std::priority_queue<event> event_queue;

	igraph_t transport_network;
	igraph_t epidemic_network;

	// PARAMATERS

	int world_size;
	int community_size;

	int community_network_degree;

	double mobility_rate;
	double community_infection_rate;
	double transport_infection_rate;
	double recovery_rate;
	double immunity_loss_rate;

	double initial_prevalence;

	bool initial_uniform_site_distribution = false;
	int initial_site = 0;

	// ---------

	std::exponential_distribution<double> exponential_distribution = std::exponential_distribution<double>(1);

	std::vector<std::vector<double>> transport_transition_matrix;
	std::vector<std::discrete_distribution<int>> transport_transition_distributions;

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

	void infer_transport_transition_matrix(void);

	void initialise_dynamics(void);

	void move_event_handler(const event &);
	void infection_event_handler(const event &);
	void recovery_event_handler(const event &);
	void immunity_loss_event_handler(const event &);

	void infection_spreading(const event &);

	std::list<std::tuple<double, int, int, int, double, double, double>> timeseries = std::list<std::tuple<double, int, int, int, double, double, double>>();

	void update_timeseries(const double &, const bool &);
	void update_timeseries(const double &);

	void print_state(void);
};
