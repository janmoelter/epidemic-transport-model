#include "epidemic_transport_model.hpp"

epidemic_transport_model_event::epidemic_transport_model_event(const double& time, const int& subject, const int& sender, const epidemic_transport_model_event_action& action)
{
	this->time = time;
	this->subject = subject;
	this->sender = sender;
	this->action = action;
}

epidemic_transport_model_event::epidemic_transport_model_event()
{
	
}

epidemic_transport_model_event::~epidemic_transport_model_event()
{
	
}

bool operator<(const epidemic_transport_model_event& e, const epidemic_transport_model_event& f) {
	return e.time > f.time;
}

std::ostream& operator<<(std::ostream& ostream, const epidemic_transport_model_event& e)
{
	ostream << std::fixed << std::setw(10) << std::setprecision(6) << e.time << " : " << e.subject;

	switch (e.action)
	{
		case epidemic_transport_model_event_action::MOVE:
			ostream << " (MOVE)";
			break;
		
		case epidemic_transport_model_event_action::INFECTION:
			ostream << " <-- " << e.sender << " (INFECTION)";
			break;

		case epidemic_transport_model_event_action::TRANSIENT_INFECTION:
			ostream << " <-- " << e.sender << " (TRANSIENT INFECTION)";
			break;

		case epidemic_transport_model_event_action::RECOVERY:
		ostream << " (RECOVERY)";
			break;

		default:
			break;
	}

	return ostream;
}


void print_attributes(const igraph_t *g) {

	igraph_strvector_t gnames, vnames, enames;
	igraph_vector_t gtypes, vtypes, etypes;

	// gnames: String vector, the names of the graph attributes.
	// gtypes: Numeric vector, the types of the graph attributes.
	// vnames: String vector, the names of the vertex attributes.
	// vtypes: Numeric vector, the types of the vertex attributes.
	// enames: String vector, the names of the edge attributes.
	// etypes: Numeric vector, the types of the edge attributes. 
	

	igraph_vector_init(&gtypes, 0);
	igraph_vector_init(&vtypes, 0);
	igraph_vector_init(&etypes, 0);
	igraph_strvector_init(&gnames, 0);
	igraph_strvector_init(&vnames, 0);
	igraph_strvector_init(&enames, 0);

	igraph_cattribute_list(g, &gnames, &gtypes, &vnames, &vtypes, &enames, &etypes);


	std::string _attrib_name;

	std::cout << "********************************************************************************" << std::endl;

	for (long int a = 0; a < igraph_strvector_size(&gnames); a++)
	{
		_attrib_name = STR(gnames, a);

		std::cout << " - " << _attrib_name << ": ";
		if (VECTOR(gtypes)[a] == IGRAPH_ATTRIBUTE_NUMERIC)
		{
			std::cout << GAN(g, _attrib_name.c_str());
		}
		else
		{
			std::cout << GAS(g, _attrib_name.c_str());
		}
	}
	std::cout << std::endl;

	std::cout << "********************************************************************************" << std::endl;

	for (long int i = 0; i < igraph_vcount(g); i++)
	{
		std::cout << "Vertex " << i << std::endl;
		for (long int a = 0; a < igraph_strvector_size(&vnames); a++) {
			_attrib_name = STR(vnames, a);

			std::cout << " - " << _attrib_name << ": ";
			if (VECTOR(vtypes)[a] == IGRAPH_ATTRIBUTE_NUMERIC)
			{
				std::cout << VAN(g, _attrib_name.c_str(), i);
			}
			else
			{
				std::cout << VAS(g, _attrib_name.c_str(), i);
			}
			std::cout << std::endl;
		}
	}

	std::cout << "********************************************************************************" << std::endl;

	for (long int i = 0; i < igraph_ecount(g); i++)
	{
		std::cout << "Edge " << (int)IGRAPH_FROM(g, i) << "--" << (int)IGRAPH_TO(g, i) << std::endl;
		for (long int a = 0; a < igraph_strvector_size(&enames); a++) {
			_attrib_name = STR(enames, a);

			std::cout << " - " << _attrib_name << ": ";
			if (VECTOR(etypes)[a] == IGRAPH_ATTRIBUTE_NUMERIC)
			{
				std::cout << EAN(g, _attrib_name.c_str(), i);
			}
			else
			{
				std::cout << EAS(g, _attrib_name.c_str(), i);
			}
			std::cout << std::endl;
		}
	}

	igraph_strvector_destroy(&enames);
	igraph_strvector_destroy(&vnames);
	igraph_strvector_destroy(&gnames);
	igraph_vector_destroy(&etypes);
	igraph_vector_destroy(&vtypes);
	igraph_vector_destroy(&gtypes);

}

epidemic_transport_model::epidemic_transport_model(const std::string& transport_network_file, const int& community_size, const int& k, const double& prevalence, const double& mobility_rate, const double& infection_rate, const double& recovery_rate, const double& coupling)
{
	igraph_set_attribute_table(&igraph_cattribute_table);

	std::mt19937 gen(this->random_number_engine());

	this->mobility_rate = mobility_rate;
	this->infection_rate = infection_rate;
	this->recovery_rate = recovery_rate;
	this->coupling = coupling;

	this->T_transport_transition = std::exponential_distribution<double>(this->mobility_rate);
	this->T_infection = std::exponential_distribution<double>(this->infection_rate);
	this->T_transient_infection = std::exponential_distribution<double>(this->coupling * this->infection_rate);
	this->T_recovery = std::exponential_distribution<double>(this->recovery_rate);


	this->community_size = community_size;
	this->prevalence = prevalence;


	// LOAD EXISTING TRANSPORT NETWORK

	FILE *ifile;

	ifile = fopen(transport_network_file.c_str(), "r");
	if (ifile == 0) {
		std::cout << "Cannot open file." << std::endl;
	}

	igraph_read_graph_graphml(&this->transport_network, ifile, 0);
	fclose(ifile);

	if (igraph_is_directed(&this->transport_network))
	{
		//igraph_to_directed(&this->transport_network, IGRAPH_TO_DIRECTED_MUTUAL);
	}

	this->world_size = igraph_vcount(&this->transport_network);

	// CREATE EPIDEMIC NETWORK (BASE LAYER) AND INITIALISE STATES

	//igraph_full(&this->epidemic_network, this->community_size, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
	igraph_k_regular_game(&this->epidemic_network, this->community_size, k, false, false);



	this->future__T_move = std::vector<double>(this->community_size, 0.);
	this->future__T_recovery = std::vector<double>(this->community_size, 0.);

	// INFER TRANSPORT NETWORK TRANSITION PROBABILITIES

	this->infer_transport_transition_matrix();

	// +++

	this->community_adjacents = std::vector<std::set<int>>(this->community_size);
	
	igraph_vs_t vs;
	igraph_vit_t vit;

	for (int n = 0; n < this->community_size; n++)
	{
		igraph_vs_adj(&vs, n, IGRAPH_OUT);
		igraph_vit_create(&this->epidemic_network, vs, &vit);

		while (!IGRAPH_VIT_END(vit)) {
			this->community_adjacents[n].insert(IGRAPH_VIT_GET(vit));
			IGRAPH_VIT_NEXT(vit);
		}
	}

	igraph_vit_destroy(&vit);
	igraph_vs_destroy(&vs);


	this->site_occupancy = std::vector<std::set<int>>(igraph_vcount(&this->transport_network));

}

epidemic_transport_model::~epidemic_transport_model()
{
	igraph_cattribute_remove_all(&this->transport_network, true, true, true);
	igraph_cattribute_remove_all(&this->epidemic_network, true, true, true);

	igraph_destroy(&this->transport_network);
	igraph_destroy(&this->epidemic_network);
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

void epidemic_transport_model::initialise_state()
{
	igraph_vs_t vs;
	igraph_vit_t vit;

	igraph_vs_all(&vs);
	igraph_vit_create(&this->epidemic_network, vs, &vit);

	while (!IGRAPH_VIT_END(vit)) {
		SETVAS(&this->epidemic_network, "health", IGRAPH_VIT_GET(vit), "S");
		SETVAN(&this->epidemic_network, "site",  IGRAPH_VIT_GET(vit), 0);

		IGRAPH_VIT_NEXT(vit);
	}

	igraph_vit_create(&this->epidemic_network, vs, &vit);
	int n;

	int x;

	double T;

	while (!IGRAPH_VIT_END(vit))
	{
		n = IGRAPH_VIT_GET(vit);

		x = std::uniform_int_distribution<int>(0, igraph_vcount(&this->transport_network)-1)(this->random_number_engine);
		//x = 0;
		SETVAN(&this->epidemic_network, "site",  n, x);
		this->site_occupancy[x].insert(n);

		T = 0.;

		if (std::uniform_real_distribution<double>(0, 1)(this->random_number_engine) < this->prevalence)
		{
			this->event_queue.emplace(T, n, NULL, epidemic_transport_model_event_action::INFECTION);
		}

		double T = 0 + this->T_transport_transition(this->random_number_engine);
		this->future__T_move[n] = T;
		if (T < this->runtime)
		{
			this->event_queue.emplace(T, n, NULL, epidemic_transport_model_event_action::MOVE);
		}

		IGRAPH_VIT_NEXT(vit);
	}

	igraph_vit_destroy(&vit);
	igraph_vs_destroy(&vs);


	this->event_queue.emplace(this->time_quantum, NULL, NULL, epidemic_transport_model_event_action::UNSPECIFIED);
}

void epidemic_transport_model::move_event_handler(const epidemic_transport_model_event& e)
{
	int x = VAN(&this->epidemic_network, "site", e.subject);
	this->site_occupancy[x].erase(e.subject);

	x = this->transport_transition_distributions[x](this->random_number_engine);

	SETVAN(&this->epidemic_network, "site", e.subject, x);
	this->site_occupancy[x].insert(e.subject);


	double T = e.time + this->T_transport_transition(this->random_number_engine);

	this->future__T_move[e.subject] = T;
	if (T < this->runtime)
	{
		this->event_queue.emplace(T, e.subject, NULL, epidemic_transport_model_event_action::MOVE);
	}

	for (auto it = this->site_occupancy[x].begin(); it != this->site_occupancy[x].end(); ++it)
	{
		if (e.subject != *it)
		{
			if (strcmp(VAS(&this->epidemic_network, "health", *it), "I") == 0)
			{
				this->infection_spreading(epidemic_transport_model_event(e.time, e.subject, *it, epidemic_transport_model_event_action::TRANSIENT_INFECTION));
			}

			if (strcmp(VAS(&this->epidemic_network, "health", e.subject), "I") == 0)
			{
				this->infection_spreading(epidemic_transport_model_event(e.time, *it, e.subject, epidemic_transport_model_event_action::TRANSIENT_INFECTION));
			}
		}
	}
}

void epidemic_transport_model::infection_event_handler(const epidemic_transport_model_event& e)
{
	if (strcmp(VAS(&this->epidemic_network, "health", e.subject), "S") == 0)
	{
		SETVAS(&this->epidemic_network, "health", e.subject, "I");

		double T = e.time + this->T_recovery(this->random_number_engine);

		this->future__T_recovery[e.subject] = T;
		if (T < this->runtime)
		{
			this->event_queue.emplace(T, e.subject, NULL, epidemic_transport_model_event_action::RECOVERY);
		}

		for (auto it = this->community_adjacents[e.subject].begin(); it != this->community_adjacents[e.subject].end(); ++it)
		{
			this->infection_spreading(epidemic_transport_model_event(e.time, *it, e.subject, epidemic_transport_model_event_action::INFECTION));
		}

		if (this->coupling > 0)
		{
			int x = VAN(&this->epidemic_network, "site", e.subject);

			for (auto it = this->site_occupancy[x].begin(); it != this->site_occupancy[x].end(); ++it)
			{
				if (e.subject != *it)
				{
					this->infection_spreading(epidemic_transport_model_event(e.time, *it, e.subject, epidemic_transport_model_event_action::TRANSIENT_INFECTION));
				}
			}
		}
	}
	else
	{
		this->infection_spreading(e);
	}
}

void epidemic_transport_model::infection_spreading(const epidemic_transport_model_event& e)
{
	switch (e.action)
	{
		case epidemic_transport_model_event_action::INFECTION:
			if (this->future__T_recovery[e.subject] < this->future__T_recovery[e.sender])
			{
				double T = std::max(e.time, this->future__T_recovery[e.subject]) + this->T_infection(this->random_number_engine);

				if (T < this->future__T_recovery[e.sender])
				{
					if (T < this->runtime)
					{
						this->event_queue.emplace(T, e.subject, e.sender, e.action);
					}
				}
			}

			break;

		case epidemic_transport_model_event_action::TRANSIENT_INFECTION:
			if (this->future__T_recovery[e.subject] < this->future__T_recovery[e.sender])
			{
				double T = std::max(e.time, this->future__T_recovery[e.subject]) + this->T_transient_infection(this->random_number_engine);

				if (T < this->future__T_recovery[e.sender] && T < std::min(this->future__T_move[e.subject], this->future__T_move[e.sender]))
				{
					//std::cout << "current time: " << e.time << std::endl;
					//std::cout << "current sites: " << (int)VAN(&this->epidemic_network, "site", e.subject) << " & " << (int)VAN(&this->epidemic_network, "site", e.sender) << std::endl;
					//std::cout << "infection time: " << T << std::endl;
					//std::cout << "moving times: " << this->future__T_move[e.subject] << " & " << this->future__T_move[e.sender] << std::endl;

					if (T < this->runtime)
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

void epidemic_transport_model::recovery_event_handler(const epidemic_transport_model_event& e)
{
	SETVAS(&this->epidemic_network, "health", e.subject, "S");
	this->future__T_recovery[e.subject] = -1.;
}


void epidemic_transport_model::update_timeseries(const double& t, const bool& std_out)
{
	int I = 0;

	igraph_vs_t vs;
	igraph_vit_t vit;

	igraph_vs_all(&vs);
	igraph_vit_create(&this->epidemic_network, vs, &vit);

	while (!IGRAPH_VIT_END(vit)) {
		if (strcmp(VAS(&this->epidemic_network, "health", IGRAPH_VIT_GET(vit)), "I") == 0)
		{
			I += 1;
		}

		IGRAPH_VIT_NEXT(vit);
	}

	igraph_vit_destroy(&vit);
	igraph_vs_destroy(&vs);


	this->timeseries.emplace_back(t, I, (double)I / this->community_size);
	
	if (std_out)
	{
		std::cout << std::fixed << std::setw(10) << std::setprecision(6) << t << "   " << std::fixed << std::setw(10) << std::setprecision(3) << (double) I / this->community_size << std::endl;
	}
}

void epidemic_transport_model::update_timeseries(const double& t)
{
	this->update_timeseries(t, false);
}

std::list<std::tuple<double, int, double>> epidemic_transport_model::run(const double& time)
{
	this->runtime = time;
	this->initialise_state();


	epidemic_transport_model_event e;

	while(!this->event_queue.empty())
	{
		e = this->event_queue.top();
		//std::cout << e << std::endl;
		
		this->event_queue.pop();

		switch (e.action)
		{
			case epidemic_transport_model_event_action::MOVE:
				this->move_event_handler(e);
				break;
			
			case epidemic_transport_model_event_action::INFECTION:
				this->infection_event_handler(e);
				break;

			case epidemic_transport_model_event_action::TRANSIENT_INFECTION:
				this->infection_event_handler(e);
				break;

			case epidemic_transport_model_event_action::RECOVERY:
				this->recovery_event_handler(e);
				break;

			default:
				break;
		}

		if (e.time > 0)
		{
			this->update_timeseries(e.time, false);
			//this->print_state();
		}
	}

	this->update_timeseries(this->runtime, false);

	return this->timeseries;
}

void epidemic_transport_model::print_state()
{
	igraph_vs_t vs;
	igraph_vit_t vit;

	igraph_vs_all(&vs);
	igraph_vit_create(&this->epidemic_network, vs, &vit);

	std::cout << "{";

	while (!IGRAPH_VIT_END(vit))
	{
		if (strcmp(VAS(&this->epidemic_network, "health", IGRAPH_VIT_GET(vit)), "S") == 0)
		{
			std::cout << " ";
		}
		else
		{
			std::cout << "*";
		}

		IGRAPH_VIT_NEXT(vit);
	}

	std::cout << "}" << std::endl;

	igraph_vit_destroy(&vit);
	igraph_vs_destroy(&vs);
}
