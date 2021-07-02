#include "epidemic_transport_model.hpp"



epidemic_transport_model::event::event(const double& time, const int& subject, const int& sender, const epidemic_transport_model::event::ACTION& action)
{
	this->time = time;
	this->subject = subject;
	this->sender = sender;
	this->action = action;
}

epidemic_transport_model::event::event()
{
	
}

epidemic_transport_model::event::~event()
{
	
}


bool operator<(epidemic_transport_model::event const& e, epidemic_transport_model::event const& f) {
	return e.time > f.time;
}

std::ostream& operator<<(std::ostream& ostream, const epidemic_transport_model::event& e)
{
	ostream << std::fixed << std::setw(10) << std::setprecision(6) << e.time << " : " << e.subject;

	switch (e.action)
	{
		case epidemic_transport_model::event::ACTION::MOVE:
			ostream << " (MOVE)";
			break;
		
		case epidemic_transport_model::event::ACTION::COMMUNITY_INFECTION:
			ostream << " <-- " << e.sender << " (INFECTION)";
			break;

		case epidemic_transport_model::event::ACTION::TRANSPORT_INFECTION:
			ostream << " <-- " << e.sender << " (TRANSIENT INFECTION)";
			break;

		case epidemic_transport_model::event::ACTION::RECOVERY:
		ostream << " (RECOVERY)";
			break;

		default:
			break;
	}

	return ostream;
}



epidemic_transport_model::epidemic_transport_model(const std::string& transport_network_file, const int& community_size, const int& community_network_degree, const double& initial_prevalence, const double& mobility_rate, const double& community_infection_rate, const double& transport_infection_rate, const double& recovery_rate)
{
	igraph_set_attribute_table(&igraph_cattribute_table);

	std::mt19937 gen(this->random_number_engine());


	// LOAD TRANSPORT NETWORK FROM FILE

	FILE *ifile;

	ifile = fopen(transport_network_file.c_str(), "r");
	if (ifile == 0) {
		std::cout << "Cannot open file." << std::endl;
	}

	igraph_read_graph_graphml(&this->transport_network, ifile, 0);
	fclose(ifile);

	//if (igraph_is_directed(&this->transport_network))
	//{
	//	igraph_to_directed(&this->transport_network, IGRAPH_TO_DIRECTED_MUTUAL);
	//}

	this->world_size = igraph_vcount(&this->transport_network);


	// SET REMAINING PARAMETERS

	this->community_size = community_size;
	this->community_network_degree = community_network_degree;
	this->initial_prevalence = initial_prevalence;

	this->mobility_rate = mobility_rate;
	this->community_infection_rate = community_infection_rate;
	this->transport_infection_rate = transport_infection_rate;
	this->recovery_rate = recovery_rate;


	// CREATE A RANDOM EPIDEMIC COMMUNITY NETWORK

	igraph_k_regular_game(&this->epidemic_network, this->community_size, this->community_network_degree, false, false);

	this->community_contacts = std::vector<std::set<int>>(this->community_size);
	
	igraph_vs_t vs;
	igraph_vit_t vit;

	for (int n = 0; n < this->community_size; n++)
	{
		igraph_vs_adj(&vs, n, IGRAPH_OUT);
		igraph_vit_create(&this->epidemic_network, vs, &vit);

		while (!IGRAPH_VIT_END(vit)) {
			this->community_contacts[n].insert(IGRAPH_VIT_GET(vit));
			IGRAPH_VIT_NEXT(vit);
		}
	}

	igraph_vit_destroy(&vit);
	igraph_vs_destroy(&vs);


	// PRE-INITIALISE STATE

	this->state_health = std::vector<epidemic_transport_model::HEALTH>(this->community_size, epidemic_transport_model::HEALTH::S);
	this->state_site = std::vector<int>(this->community_size, 0);

	this->site_occupancy = std::vector<std::set<int>>(this->world_size, std::set<int>());
	for (int n = 0; n < this->community_size; n++)
	{
		this->site_occupancy[0].insert(n);
	}
	
	this->state_site_future_transition_time = std::vector<double>(this->community_size, 0.);
	this->state_health_future_recovery_time = std::vector<double>(this->community_size, 0.);
}

epidemic_transport_model::~epidemic_transport_model()
{
	igraph_cattribute_remove_all(&this->transport_network, true, true, true);
	igraph_cattribute_remove_all(&this->epidemic_network, true, true, true);

	igraph_destroy(&this->transport_network);
	igraph_destroy(&this->epidemic_network);
}


std::list<std::tuple<double, int, double>> epidemic_transport_model::simulate(const double& time, const bool& verbose)
{
	this->simulation_time = time;
	this->initialise_dynamics();


	event e;

	while(!this->event_queue.empty())
	{
		e = this->event_queue.top();
		this->event_queue.pop();

		switch (e.action)
		{
			case epidemic_transport_model::event::ACTION::MOVE:
				this->move_event_handler(e);
				break;
			
			case epidemic_transport_model::event::ACTION::COMMUNITY_INFECTION:
			case epidemic_transport_model::event::ACTION::TRANSPORT_INFECTION:
				this->infection_event_handler(e);
				break;

			case epidemic_transport_model::event::ACTION::RECOVERY:
				this->recovery_event_handler(e);
				break;

			default:
				break;
		}

		if (e.time > 0)
		{
			this->update_timeseries(e.time, false);
			if (verbose)
			{
				this->print_state();
			}
		}
	}

	this->update_timeseries(this->simulation_time, false);

	return this->timeseries;
}

std::list<std::tuple<double, int, double>> epidemic_transport_model::simulate(const double& time)
{
	return epidemic_transport_model::simulate(time, false);
}


void epidemic_transport_model::infer_transport_transition_matrix()
{
	this->transport_transition_matrix = std::vector<std::vector<double>>(this->world_size, std::vector<double>(this->world_size, 0.));

	this->transport_transition_distributions = std::vector<std::discrete_distribution<int>>(this->world_size);

	bool is_weighted = igraph_cattribute_has_attr(&this->transport_network, IGRAPH_ATTRIBUTE_EDGE, "weight");

	std::tuple<int, int> arc;
	double weight;

	for (long int a = 0; a < igraph_ecount(&this->transport_network); a++)
	{
		arc = std::make_tuple((int)IGRAPH_FROM(&this->transport_network, a), (int)IGRAPH_TO(&this->transport_network, a));

		this->transport_transition_matrix[std::get<0>(arc)][std::get<1>(arc)] = 1.;
		this->transport_transition_matrix[std::get<1>(arc)][std::get<0>(arc)] = 1.;
		if (is_weighted)
		{
			weight = EAN(&this->transport_network, "weight", a);
			if(weight < 0)
			{
				std::cout << "Transport network shows negative weights. These will be set to 0.";
				weight = 0;
			}

			this->transport_transition_matrix[std::get<0>(arc)][std::get<1>(arc)] *= weight;
			this->transport_transition_matrix[std::get<1>(arc)][std::get<0>(arc)] *= weight;
		}
		
	}

	double weight_sum;

	for (int x = 0; x < this->world_size; x++)
	{
		weight_sum = std::accumulate(this->transport_transition_matrix[x].begin(), this->transport_transition_matrix[x].end(), 0.);

		if (weight_sum > 0.)
		{
			std::for_each(this->transport_transition_matrix[x].begin(), this->transport_transition_matrix[x].end(), [weight_sum](double &w){ w /= weight_sum; });
		}

		//for (int i = 0; i < this->transport_transition_matrix[x].size(); i++)
		//{
		//	std::cout << std::fixed << std::setw(7) << std::setprecision(3) << this->transport_transition_matrix[x][i];
		//}
		//std::cout << std::endl;

		this->transport_transition_distributions[x] = std::discrete_distribution<int>(this->transport_transition_matrix[x].begin(), this->transport_transition_matrix[x].end());
	}
}

void epidemic_transport_model::initialise_dynamics()
{

	this->infer_transport_transition_matrix();


	int x;
	epidemic_transport_model::HEALTH h;

	double T;


	std::uniform_int_distribution<int> initial_site_distribution;

	if (this->initial_uniform_site_distribution)
	{
		initial_site_distribution = std::uniform_int_distribution<int>(0, this->world_size-1);
	}
	else
	{
		initial_site_distribution = std::uniform_int_distribution<int>(this->initial_site, this->initial_site);
	}

	for (int n = 0; n < this->community_size; n++)
	{

		x = 0;
		this->site_occupancy[x].erase(n);
		x = initial_site_distribution(this->random_number_engine);
		this->state_site[n] = x;
		this->site_occupancy[x].insert(n);


		T = 0.;

		if (std::uniform_real_distribution<double>(0, 1)(this->random_number_engine) < this->initial_prevalence)
		{
			this->event_queue.emplace(T, n, NULL, epidemic_transport_model::event::ACTION::COMMUNITY_INFECTION);
		}

		T = 0 + this->exponential_distribution(this->random_number_engine) / this->mobility_rate;
		this->state_site_future_transition_time[n] = T;
		if (T < this->simulation_time)
		{
			this->event_queue.emplace(T, n, NULL, epidemic_transport_model::event::ACTION::MOVE);
		}

	}

	T = std::numeric_limits<double>::min();
	this->event_queue.emplace(T, NULL, NULL, epidemic_transport_model::event::ACTION::UNSPECIFIED);
}

void epidemic_transport_model::move_event_handler(const event& e)
{
	int x;

	x = this->state_site[e.subject];
	this->site_occupancy[x].erase(e.subject);
	x = this->transport_transition_distributions[x](this->random_number_engine);
	this->state_site[e.subject] = x;
	this->site_occupancy[x].insert(e.subject);


	double T = e.time + this->exponential_distribution(this->random_number_engine) / this->mobility_rate;
	this->state_site_future_transition_time[e.subject] = T;
	if (T < this->simulation_time)
	{
		this->event_queue.emplace(T, e.subject, NULL, epidemic_transport_model::event::ACTION::MOVE);
	}


	for (auto& a : this->site_occupancy[x])
	{
		if (e.subject != a)
		{
			if (this->state_health[e.subject] == epidemic_transport_model::HEALTH::I)
			{
				this->infection_spreading(event(e.time, a, e.subject, epidemic_transport_model::event::ACTION::TRANSPORT_INFECTION));
			}

			if (this->state_health[a] == epidemic_transport_model::HEALTH::I)
			{
				this->infection_spreading(event(e.time, e.subject, a, epidemic_transport_model::event::ACTION::TRANSPORT_INFECTION));
			}
		}
	}
}

void epidemic_transport_model::infection_event_handler(const event& e)
{
	double T;

	if (this->state_health[e.subject] == epidemic_transport_model::HEALTH::S)
	{
		this->state_health[e.subject] = epidemic_transport_model::HEALTH::I;


		T = e.time + this->exponential_distribution(this->random_number_engine) / this->recovery_rate;
		this->state_health_future_recovery_time[e.subject] = T;
		if (T < this->simulation_time)
		{
			this->event_queue.emplace(T, e.subject, NULL, epidemic_transport_model::event::ACTION::RECOVERY);
		}

		for (auto& n : this->community_contacts[e.subject])
		{
			this->infection_spreading(event(e.time, n, e.subject, epidemic_transport_model::event::ACTION::COMMUNITY_INFECTION));
		}

		for (auto& n : this->site_occupancy[this->state_site[e.subject]])
		{
			if (e.subject != n)
			{
				this->infection_spreading(event(e.time, n, e.subject, epidemic_transport_model::event::ACTION::TRANSPORT_INFECTION));
			}
		}
	}
	else
	{
		this->infection_spreading(e);
	}
}

void epidemic_transport_model::infection_spreading(const event& e)
{
	double T;

	double infection_rate;
	double T_infection_path_breaking;

	switch (e.action)
	{
		case epidemic_transport_model::event::ACTION::COMMUNITY_INFECTION:
			infection_rate = this->community_infection_rate;
			T_infection_path_breaking = std::numeric_limits<double>::infinity();

			break;

		case epidemic_transport_model::event::ACTION::TRANSPORT_INFECTION:
			infection_rate = this->transport_infection_rate;
			T_infection_path_breaking = std::min(this->state_site_future_transition_time[e.subject], this->state_site_future_transition_time[e.sender]);

			break;

		default:
			break;
	}


	switch (e.action)
	{
		case epidemic_transport_model::event::ACTION::COMMUNITY_INFECTION:
		case epidemic_transport_model::event::ACTION::TRANSPORT_INFECTION:

			if (this->state_health_future_recovery_time[e.subject] < this->state_health_future_recovery_time[e.sender])
			{
				T = std::max(e.time, this->state_health_future_recovery_time[e.subject]) + this->exponential_distribution(this->random_number_engine) / infection_rate;
				
				if (T < this->state_health_future_recovery_time[e.sender] && T < T_infection_path_breaking)
				{
					if (T < this->simulation_time)
					{
						this->event_queue.emplace(T, e.subject, e.sender, e.action);
					}
				}
			}

			break;

		default:
			break;
	}
}

void epidemic_transport_model::recovery_event_handler(const event& e)
{
	this->state_health[e.subject] = epidemic_transport_model::HEALTH::S;
	this->state_health_future_recovery_time[e.subject] = -std::numeric_limits<double>::infinity();
}


void epidemic_transport_model::update_timeseries(const double& t, const bool& std_out)
{
	int I = std::count(this->state_health.begin(), this->state_health.end(), epidemic_transport_model::HEALTH::I);

	this->timeseries.emplace_back(t, I, (double)I / this->community_size);
	
	if (std_out)
	{
		std::cout << std::fixed << std::setw(10) << std::setprecision(6) << t << "   " << std::fixed << std::setw(10) << std::setprecision(3) << (double) I / this->community_size << std::endl;
	}
}

void epidemic_transport_model::update_timeseries(const double& t)
{
	epidemic_transport_model::update_timeseries(t, false);
}

void epidemic_transport_model::print_state()
{
	std::cout << "{";

	for (int n = 0; n < this->community_size; n++)
	{
		if (this->state_health[n] == epidemic_transport_model::HEALTH::S)
		{
			std::cout << " ";
		}
		else
		{
			std::cout << "*";
		}
	}

	std::cout << "}" << std::endl;
}
