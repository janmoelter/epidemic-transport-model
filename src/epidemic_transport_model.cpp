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
	if (e.time != f.time)
	{
		return e.time > f.time;
	}
	else
	{
		return (e.time >= f.time) && (e.action >= f.action);
	}
}

std::ostream& operator<<(std::ostream& ostream, const epidemic_transport_model::event& e)
{
	ostream << std::fixed << "[" << std::setw(10) << std::setprecision(6) << e.time << "]";

	switch (e.action)
	{
		case epidemic_transport_model::event::ACTION::MOVE:
			ostream << " MOVE : " << e.subject;
			break;
		
		case epidemic_transport_model::event::ACTION::COMMUNITY_INFECTION:
			ostream << " COMMUNITY INFECTION: " << e.subject << " <--- " << e.sender;
			break;

		case epidemic_transport_model::event::ACTION::TRANSPORT_INFECTION:
			ostream << " TRANSPORT INFECTION: " << e.subject << " <--- " << e.sender;
			break;

		case epidemic_transport_model::event::ACTION::RECOVERY:
			ostream << " RECOVERY : " << e.subject;
			break;

		case epidemic_transport_model::event::ACTION::IMMUNITY_LOSS:
			ostream << " IMMUNITY LOSS : " << e.subject;
			break;

		default:
			break;
	}

	return ostream;
}



void epidemic_transport_model::initialize(std::array<igraph_t *,2>& transport_networks, const std::function<double(const double&)>& transport_network_interpolation_functional, const int& community_size, const int& community_network_degree, const int& initial_site, const double& initial_prevalence, const double& mobility_rate, const double& community_infection_rate, const double& transport_infection_rate, const double& recovery_rate, const double& immunity_loss_rate)
{
	std::mt19937 gen(this->random_number_engine());


	// COPY TRANSPORT NETWORKS

	for (size_t r = 0; r < transport_networks.size(); r++)
	{
		igraph_copy(&this->transport_networks[r], transport_networks[r]);
	}


	// SET REMAINING PARAMETERS

	this->world_size = igraph_vcount(&this->transport_networks[0]);
	if (igraph_vcount(&this->transport_networks[1]) != this->world_size)
	{
		throw std::invalid_argument("The number of sites in the two transport networks has to be the same.");
	}

	this->transport_network_interpolation_functional = transport_network_interpolation_functional;

	this->community_size = community_size;
	this->community_network_degree = community_network_degree;

	this->initial_site = initial_site;
	this->initial_prevalence = initial_prevalence;

	this->mobility_rate = mobility_rate;
	this->community_infection_rate = community_infection_rate;
	this->transport_infection_rate = transport_infection_rate;
	this->recovery_rate = recovery_rate;
	this->immunity_loss_rate = immunity_loss_rate;


	// CREATE A RANDOM EPIDEMIC COMMUNITY NETWORK

	igraph_k_regular_game(&this->epidemic_network, this->community_size, this->community_network_degree, false, false);

	this->community_contacts = std::vector<std::set<int>>(this->community_size);
	
	igraph_vs_t vs;
	igraph_vit_t vit;

	for (size_t n = 0; n < this->community_size; n++)
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
	for (size_t n = 0; n < this->community_size; n++)
	{
		this->site_occupancy[0].insert(n);
	}
	
	this->state_site_future_transition_time = std::vector<double>(this->community_size, -std::numeric_limits<double>::infinity());
	this->state_health_future_recovery_time = std::vector<double>(this->community_size, -std::numeric_limits<double>::infinity());
	this->state_health_future_immunity_loss_time = std::vector<double>(this->community_size, -std::numeric_limits<double>::infinity());
}

epidemic_transport_model::epidemic_transport_model(std::array<igraph_t *,2>& transport_networks, const std::function<double(const double&)>& transport_network_interpolation_functional, const int& community_size, const int& community_network_degree, const int& initial_site, const double& initial_prevalence, const double& mobility_rate, const double& community_infection_rate, const double& transport_infection_rate, const double& recovery_rate, const double& immunity_loss_rate)
{
	this->initialize(transport_networks, transport_network_interpolation_functional, community_size, community_network_degree, initial_site, initial_prevalence, mobility_rate, community_infection_rate, transport_infection_rate, recovery_rate, immunity_loss_rate);
}

epidemic_transport_model::epidemic_transport_model(igraph_t* transport_network, const int& community_size, const int& community_network_degree, const int& initial_site, const double& initial_prevalence, const double& mobility_rate, const double& community_infection_rate, const double& transport_infection_rate, const double& recovery_rate, const double& immunity_loss_rate)
{
	std::array<igraph_t *,2> transport_networks = {transport_network, transport_network};
	std::function<double(const double&)> transport_network_interpolation_functional = [] (double t)->double { return 0.; };

	this->initialize(transport_networks, transport_network_interpolation_functional, community_size, community_network_degree, initial_site, initial_prevalence, mobility_rate, community_infection_rate, transport_infection_rate, recovery_rate, immunity_loss_rate);
}

epidemic_transport_model::epidemic_transport_model()
{

}

epidemic_transport_model::~epidemic_transport_model()
{
	for (size_t r = 0; r < this->transport_networks.size(); r++)
	{
		igraph_cattribute_remove_all(&this->transport_networks[r], true, true, true);
		igraph_destroy(&this->transport_networks[r]);
	}

	igraph_cattribute_remove_all(&this->epidemic_network, true, true, true);
	igraph_destroy(&this->epidemic_network);
}


std::list<std::tuple<double, int, int, int, double, double, double>> epidemic_transport_model::simulate(const double& time, const bool& verbose)
{
	this->simulation_time = time;
	this->initialise_dynamics();


	event e;

	while(!this->event_queue.empty())
	{
		e = this->event_queue.top();
		this->event_queue.pop();

		if (verbose)
		{
			std::cout << e << std::endl;
		}

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

			case epidemic_transport_model::event::ACTION::IMMUNITY_LOSS:
				this->immunity_loss_event_handler(e);
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

std::list<std::tuple<double, int, int, int, double, double, double>> epidemic_transport_model::simulate(const double& time)
{
	return epidemic_transport_model::simulate(time, false);
}


std::ostream& operator<<(std::ostream& ostream, const epidemic_transport_model& _)
{
	//ostream << std::fixed << "[" << std::setw(10) << std::setprecision(6) << e.time << "]";

	ostream << std::setfill('*') << std::setw(80) << "" << std::endl;
	ostream << "* transport network : " << "*" << std::endl;
	ostream << "* community size : " << _.community_size << std::endl;
	ostream << "* community degree : " << _.community_network_degree << std::endl;
	ostream << "* initial prevalence : " << _.initial_prevalence << std::endl;
	ostream << "* mobility rate : " << _.mobility_rate << std::endl;
	ostream << "* community infection rate : " << _.community_infection_rate << std::endl;
	ostream << "* transport infection rate : " << _.transport_infection_rate << std::endl;
	ostream << "* recovery rate : " << _.recovery_rate << std::endl;
	ostream << "* immunity loss rate : " << _.immunity_loss_rate << std::endl;
	ostream << std::setfill('*') << std::setw(80) << "" << std::endl;
	ostream << std::setfill(' ');

	if (_.timeseries.size() > 0)
	{
		int w_T = (int)ceil(log10(_.simulation_time));
		int w_N = (int)ceil(log10(_.community_size));

		ostream << std::fixed << std::right << std::setw(w_T + 1 + 9) << "time" << '\t';

		ostream << std::setw(w_N) << "#S" << '\t';
		ostream << std::setw(w_N) << "#I" << '\t';
		ostream << std::setw(w_N) << "#R" << '\t';
			
		ostream << std::setw(1 + 1 + 9) << "S" << '\t';
		ostream << std::setw(1 + 1 + 9) << "I" << '\t';
		ostream << std::setw(1 + 1 + 9) << "R" << std::endl;

		for (auto& e: _.timeseries)
		{
			ostream << std::fixed << std::right << std::setw(w_T + 1 + 9) << std::setprecision(9) << std::get<0>(e) << '\t';

			ostream << std::setw(w_N) << std::get<1>(e) << '\t';
			ostream << std::setw(w_N) << std::get<2>(e) << '\t';
			ostream << std::setw(w_N) << std::get<3>(e) << '\t';
			
			ostream << std::setw(1 + 1 + 9) << std::setprecision(9) << std::get<4>(e) << '\t';
			ostream << std::setw(1 + 1 + 9) << std::setprecision(9) << std::get<5>(e) << '\t';
			ostream << std::setw(1 + 1 + 9) << std::setprecision(9) << std::get<6>(e) << std::endl;
		}

	}

	return ostream;
}


void epidemic_transport_model::infer_transport_transition_distributions()
{
	std::vector<std::vector<double>> transition_matrix;
	bool is_directed;
	bool is_weighted;

	std::tuple<int, int> arc;
	double weight;

	for (size_t r = 0; r < this->transport_networks.size(); r++)
	{
		transition_matrix = std::vector<std::vector<double>>(this->world_size, std::vector<double>(this->world_size, 0.));
		this->transport_transition_distributions[r] = std::vector<std::discrete_distribution<int>>(this->world_size);

		is_directed = igraph_is_directed(&this->transport_networks[r]);
		is_weighted = igraph_cattribute_has_attr(&this->transport_networks[r], IGRAPH_ATTRIBUTE_EDGE, "weight");


		for (long int a = 0; a < igraph_ecount(&this->transport_networks[r]); a++)
		{
			arc = std::make_tuple((int)IGRAPH_FROM(&this->transport_networks[r], a), (int)IGRAPH_TO(&this->transport_networks[r], a));


			if (is_weighted)
			{
				weight = EAN(&this->transport_networks[r], "weight", a);
				if(weight < 0)
				{
					std::cout << "Transport network has negative weights. These will be set to 0.";
					weight = 0.;
				}

				transition_matrix[std::get<0>(arc)][std::get<1>(arc)] = weight;
			}
			else
			{
				transition_matrix[std::get<0>(arc)][std::get<1>(arc)] = 1.;
			}

			if (!is_directed)
			{
				transition_matrix[std::get<1>(arc)][std::get<0>(arc)] = transition_matrix[std::get<0>(arc)][std::get<1>(arc)];
			}
		}

		for (size_t x = 0; x < this->world_size; x++)
		{
			this->transport_transition_distributions[r][x] = std::discrete_distribution<int>(transition_matrix[x].begin(), transition_matrix[x].end());

			//for (size_t y = 0; y < this->world_size; y++)
			//{
			//	std::cout << std::fixed << std::setw(6) << std::setprecision(2) << transition_matrix[x][y];
			//}
			//std::cout << std::endl;
			
		}
	}


	//std::vector<double> P;
	//for (size_t x = 0; x < this->world_size; x++)
	//{
	//	for (size_t r = 0; r < 2; r++)
	//	{
	//		P = this->transport_transition_distributions[r][x].probabilities();
    //		for(auto p : P)
    //    		std::cout << std::fixed << std::setw(6) << std::setprecision(3) << p;
    //		std::cout << std::endl;
	//	}
	//	std::cout << "****************************************************************************************************" << std::endl;
	//}
}

void epidemic_transport_model::initialise_dynamics()
{

	this->infer_transport_transition_distributions();


	int x;
	epidemic_transport_model::HEALTH h;

	double T;


	std::uniform_int_distribution<int> initial_site_distribution;

	if (0 <= this->initial_site && this->initial_site < this->world_size)
	{
		initial_site_distribution = std::uniform_int_distribution<int>(this->initial_site, this->initial_site);
	}
	else
	{
		initial_site_distribution = std::uniform_int_distribution<int>(0, this->world_size-1);
	}

	for (size_t n = 0; n < this->community_size; n++)
	{

		x = this->state_site[n];
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
	size_t r = (this->uniform_distribution(this->random_number_engine) > this->transport_network_interpolation_functional(e.time)) ? 0 : 1;
	size_t x;

	x = this->state_site[e.subject];
	this->site_occupancy[x].erase(e.subject);
	x = this->transport_transition_distributions[r][x](this->random_number_engine);
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


		// DETERMINE FUTURE RECOVERY (I --> R) TIME
		T = e.time + this->exponential_distribution(this->random_number_engine) / this->recovery_rate;
		this->state_health_future_recovery_time[e.subject] = T;
		if (T < this->simulation_time)
		{
			this->event_queue.emplace(T, e.subject, NULL, epidemic_transport_model::event::ACTION::RECOVERY);
		}

		// DETERMINE FUTURE IMMUNITY LOSS (R --> S) TIME
		T = T + this->exponential_distribution(this->random_number_engine) / this->immunity_loss_rate;
		this->state_health_future_immunity_loss_time[e.subject] = T;
		if (T < this->simulation_time)
		{
			this->event_queue.emplace(T, e.subject, NULL, epidemic_transport_model::event::ACTION::IMMUNITY_LOSS);
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

			if (this->state_health_future_immunity_loss_time[e.subject] < this->state_health_future_recovery_time[e.sender])
			{
				T = std::max(e.time, this->state_health_future_immunity_loss_time[e.subject]) + this->exponential_distribution(this->random_number_engine) / infection_rate;
				
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
	this->state_health[e.subject] = epidemic_transport_model::HEALTH::R;
	this->state_health_future_recovery_time[e.subject] = -std::numeric_limits<double>::infinity();
}

void epidemic_transport_model::immunity_loss_event_handler(const event& e)
{
	this->state_health[e.subject] = epidemic_transport_model::HEALTH::S;
	this->state_health_future_immunity_loss_time[e.subject] = -std::numeric_limits<double>::infinity();
}


void epidemic_transport_model::update_timeseries(const double& t, const bool& std_out)
{
	int S = std::count(this->state_health.begin(), this->state_health.end(), epidemic_transport_model::HEALTH::S);
	int I = std::count(this->state_health.begin(), this->state_health.end(), epidemic_transport_model::HEALTH::I);
	int R = std::count(this->state_health.begin(), this->state_health.end(), epidemic_transport_model::HEALTH::R);

	double s = (double)S / this->community_size;
	double i = (double)I / this->community_size;
	double r = (double)R / this->community_size;

	this->timeseries.emplace_back(t, S, I, R, s, i, r);
	
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

	for (size_t n = 0; n < this->community_size; n++)
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
