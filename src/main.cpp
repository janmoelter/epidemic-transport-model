#include <iostream>
#include <fstream>

#include <getopt.h>
#include <map>

#include <regex>

#include <chrono>

#include "epidemic_transport_model.hpp"



void print_version()
{
	std::cout << "Epidemic-Transport-Model Simulation (Compilation: " << __DATE__ << " @ " << __TIME__ << ")" << std::endl;
}

void print_help()
{
	print_version();

	std::cout << "" << std::endl;
	std::cout << "Usage: simulaton -w FILE ..." << std::endl;
	std::cout << "   or  simulaton -w FILE -w FILE -f FUNCTION ..." << std::endl;

	std::cout << "" << std::endl;
	std::cout << "Options:" << std::endl;
	//           "****************************************************************************************************"
	std::cout << "  -w, --transport-network-file=FILE                           File containing a transport network in" << std::endl;
	std::cout << "                                                              the GraphML-format.                   " << std::endl;
	std::cout << "  -f, --transport-network-interpolation-function=DESCRIPTION  Function to interpolate between two   " << std::endl;
	std::cout << "                                                              transport networks (if provided). If  " << std::endl;
	std::cout << "                                                              empty, the constant zero-function is  " << std::endl;
	std::cout << "                                                              used. Other options are 'theta():t0'  " << std::endl;
	std::cout << "                                                              for a Heaviside-theta pulse starting  " << std::endl;
	std::cout << "                                                              at t0, as well as 'square(w):t0++T'   " << std::endl;
	std::cout << "                                                              and 'sine(w):t0++T' for periodic      " << std::endl;
	std::cout << "                                                              square- and sine-pulses with a width  " << std::endl;
	std::cout << "                                                              w, starting at time t0 and period-    " << std::endl;
	std::cout << "                                                              length T. Required only if two        " << std::endl;
	std::cout << "                                                              transport network files are provided. " << std::endl;
	std::cout << "  -N, --community-size=INTEGER                                Size of the community.                " << std::endl;
	std::cout << "  -k, --community-degree=INTEGER                              Degree of links in the community.     " << std::endl;
	std::cout << "  -m, --mobility-rate=NUMBER                                  Mobility rate of the random walk of   " << std::endl;
	std::cout << "                                                              individuals on the transport network. " << std::endl;
	std::cout << "  -B, --community-infection-rate=NUMBER                       Infection rate in the community layer " << std::endl;
	std::cout << "                                                              of the epidemic network.              " << std::endl;
	std::cout << "  -b, --transport-infection-rate=NUMBER                       Infection rate in the transport layer " << std::endl;
	std::cout << "                                                              of the epidemic network.              " << std::endl;
	std::cout << "  -g, --recovery-rate=NUMBER                                  Recovery rate on the epidemic network." << std::endl;
	std::cout << "  -s, --immunity-loss-rate=NUMBER                             Immunity loss rate on the epidemic    " << std::endl;
	std::cout << "                                                              network. For SIR- or SIS-epidemics    " << std::endl;
	std::cout << "                                                              this is '0' or 'inf', respectively.   " << std::endl;
	std::cout << "  -T, --time=NUMBER                                           Time-span of the simulation.          " << std::endl;
	std::cout << "  -x, --initial-site=INTEGER                                  Initial site of individuals in the    " << std::endl;
	std::cout << "                                                              transport network. Invalid values,    " << std::endl;
	std::cout << "                                                              result in a uniform distribution.     " << std::endl;
	std::cout << "                                                              Default is '-1'.                      " << std::endl;
	std::cout << "  -p, --initial-prevalence=NUMBER                             Initial prevalence of the epidemic in " << std::endl;
	std::cout << "                                                              the community.                        " << std::endl;
	std::cout << "  -a, --fractional-exponent=NUMBER                            Exponent of the fractional dynamics.  " << std::endl;
	std::cout << "                                                              Default is '1'.                       " << std::endl;
	std::cout << "  -o, --output-file=FILE                                      Output file for the simulation time-  " << std::endl;
	std::cout << "                                                              series. If this is empty, standard    " << std::endl;
	std::cout << "                                                              output is used. Default is ''.        " << std::endl;
	std::cout << "  -v, --verbose                                               " << std::endl;
	std::cout << "  -h, --help                                                  " << std::endl;
}

std::map<std::string, std::vector<std::string>> argparse(int argc, char** argv)
{
	std::map<std::string, std::vector<std::string>> argm { {"transport-network-file", {}}, {"transport-network-interpolation-function", {}}, {"community-size", {}}, {"community-degree", {}}, {"mobility-rate", {}}, {"community-infection-rate", {}}, {"transport-infection-rate", {}}, {"recovery-rate", {}}, {"immunity-loss-rate", {}}, {"time", {}}, {"initial-site", {"-1"}}, {"initial-prevalence", {}}, {"fractional-exponent", {"1"}}, {"output-file", {""}}, {"verbose", {"0"}}, {"help", {"0"}}, {"version", {"0"}}, };

	const char* const short_opts = "w:f:N:k:m:B:b:g:s:T:x:p:a:o:vhV";
	const option long_opts[] = {
		{"transport-network-file", required_argument, nullptr, 'w'},
		{"transport-network-interpolation-function", required_argument, nullptr, 'f'},
		{"community-size", required_argument, nullptr, 'N'},
		{"community-degree", required_argument, nullptr, 'k'},
		{"mobility-rate", required_argument, nullptr, 'm'},
		{"community-infection-rate", required_argument, nullptr, 'B'},
		{"transport-infection-rate", required_argument, nullptr, 'b'},
		{"recovery-rate", required_argument, nullptr, 'g'},
		{"immunity-loss-rate", required_argument, nullptr, 's'},
		{"time", required_argument, nullptr, 'T'},
		{"initial-site", required_argument, nullptr, 'x'},
		{"initial-prevalence", required_argument, nullptr, 'p'},
		{"fractional-exponent", required_argument, nullptr, 'a'},
		{"output-file", required_argument, nullptr, 'o'},
		{"verbose", no_argument, nullptr, 'v'},
		{"help", no_argument, nullptr, 'h'},
		{"version", no_argument, nullptr, 'V'},
		{nullptr, no_argument, nullptr, 0}
	};

	while (true)
	{
		const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

		if (opt == -1)
			break;

		switch (opt)
		{
			case 'w':
				argm["transport-network-file"].push_back(optarg);
				break;
			case 'f':
				argm["transport-network-interpolation-function"].push_back(optarg);
				break;
			case 'N':
				argm["community-size"].push_back(optarg);
				break;
			case 'k':
				argm["community-degree"].push_back(optarg);
				break;
			case 'm':
				argm["mobility-rate"].push_back(optarg);
				break;
			case 'B':
				argm["community-infection-rate"].push_back(optarg);
				break;
			case 'b':
				argm["transport-infection-rate"].push_back(optarg);
				break;
			case 'g':
				argm["recovery-rate"].push_back(optarg);
				break;
			case 's':
				argm["immunity-loss-rate"].push_back(optarg);
				break;
			case 'T':
				argm["time"].push_back(optarg);
				break;
			case 'x':
				argm["initial-site"].push_back(optarg);
				break;
			case 'p':
				argm["initial-prevalence"].push_back(optarg);
				break;
			case 'a':
				argm["fractional-exponent"].push_back(optarg);
				break;
			case 'o':
				argm["output-file"].push_back(optarg);
				break;
			case 'v':
				argm["verbose"].push_back("1");
				break;
			case 'V':
				argm["version"].push_back("1");
				break;
			case 'h':
			case '?':
				argm["help"].push_back("1");
				break;
			default:
				break;
		}
	}

	if (argm["transport-network-file"].size() < 2)
	{
		argm["transport-network-interpolation-function"].push_back("");
	}

	if (argm["transport-network-file"].size() > 2)
	{
		argm["help"].push_back("1");
	}

	for (auto& [opt, arg] : argm)
	{
		if (arg.size() == 0)
		{
			argm["help"].push_back("1");
			break;
		}
	}

	return argm;
}

std::function<double(const double&)> transport_network_interpolation_function_parsing(std::string transport_network_interpolation_function_string)
{
	std::regex function_string_regex;
	std::smatch matches;

	// Empty function string "".

	function_string_regex = std::regex(R"(^$)");
	if (std::regex_match(transport_network_interpolation_function_string, matches, function_string_regex))
	{
		return [] (double t)->double { return 0.; };
	}

	// Function string 'F():t0' for a functions where F defines the function
	// prototype and t0 its onset.
	//
	// For example:
	//
	//      
	//   1 ╴      ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 
	//            ┃                                      
	//   0 ━━━━━━━┩┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄
	//           t0
	//

	function_string_regex = std::regex(R"(^(\w*)\(\)\:([+-]?\d*(?:\.(?:[0-9]+))?)$)");
	if (std::regex_match(transport_network_interpolation_function_string, matches, function_string_regex))
	{
		std::string f_type = matches[1].str();
		double arg_offset = stod(matches[2].str());

		std::function<double(const double&)> f;

		if (!f_type.compare("theta"))
		{
			return [arg_offset] (double t)->double { return 1. * (t > arg_offset); };
		}
	}

	// Function string "F(w):t0++T" for a periodic function where F defines the
	// function prototype, w its width, t0 the first onset and T the function's
	// period length.
	//
	// For example:
	//
	//      
	//   1 ╴      ┏━━━━━━━━━━┓       ┏━━━━━━━━━━┓        
	//            ┃          ┃       ┃          ┃        
	//   0 ━━━━━━━┩┄┄┄┄┄┄┄┄┄┄┡━━━━━━━┩┄┄┄┄┄┄┄┄┄┄┡━━━━━━━━
	//           t0         t0+w
	//
	//            │<────── T ───────>│
	//

	function_string_regex = std::regex(R"(^(\w*)\(([+]?\d*(?:\.(?:[0-9]+))?)\)\:([+-]?\d*(?:\.(?:[0-9]+))?)\+\+((?:[+]?\d*(?:\.(?:[0-9]+))?)|inf)$)");
	if (std::regex_match(transport_network_interpolation_function_string, matches, function_string_regex))
	{
		std::string f_type = matches[1].str();
		double arg_offset = stod(matches[3].str());
		double f_width = stod(matches[2].str());
		double f_period = stod(matches[4].str());

		std::function<double(const double&)> f;

		if (!f_type.compare("square"))
		{
			f = [] (double t)->double { return (0. < t && t < 1); };
		}

		if (!f_type.compare("sine"))
		{
			f = [] (double t)->double { return pow(sin(M_PI * t), 2.); };
		}

		return [f, arg_offset, f_width, f_period] (double t)->double { return f(fmin(fmod(t - arg_offset, f_period) / f_width, 1.)) * (t > arg_offset); };
	}

	throw std::invalid_argument(transport_network_interpolation_function_string);
}

std::vector<int> community_network_degree_parsing(std::string community_network_degree_string)
{
	std::regex string_regex;
	std::smatch matches;

	string_regex = std::regex(R"(^(?:(?:\d{1,}\*)?\d{1,},)*(?:\d{1,}\*)?\d{1,}$)");
	if (std::regex_match(community_network_degree_string, matches, string_regex))
	{
		std::vector<int> community_network_degree;

		string_regex = std::regex(R"(^(?:(\d{1,})\*)?(\d{1,})$)");
		int m;
		int k;

		std::stringstream _stringstream(community_network_degree_string);
    	std::string item;
    	while (std::getline(_stringstream, item, ','))
		{
    	    if (item.length() > 0)
			{
				if (std::regex_match(item, matches, string_regex))
				{
					m = matches[1].length() > 0 ? stoi(matches[1].str()) : 1;
					k = stoi(matches[2].str());

					community_network_degree.insert(community_network_degree.end(), m, k);
				}
    	    }
    	}

		return community_network_degree;
	}

	throw std::invalid_argument(community_network_degree_string);
}

void igraph_read_graph_graphmlfile(igraph_t* graph, std::string filename)
{
	FILE *ifile;

	ifile = fopen(filename.c_str(), "r");
	if (ifile == 0)
	{
		std::cout << "Cannot open transport network file." << std::endl;
	}

	igraph_read_graph_graphml(graph, ifile, 0);
	fclose(ifile);
}


int main(int argc, char **argv)
{
	auto start_time = std::chrono::steady_clock::now();

	igraph_set_attribute_table(&igraph_cattribute_table);

	//std::cout << "+ " << __LINE__ << std::endl;

	std::map<std::string, std::vector<std::string>> argm = argparse(argc, argv);


	bool version = (bool)std::stoi(argm["version"].back());
	if (version == true)
	{
		print_version();
		return 0;
	}

	bool help = (bool)std::stoi(argm["help"].back());
	if (help == true)
	{
		print_help();
		return 0;
	}

	std::vector<std::string> transport_network_files = argm["transport-network-file"];
	std::function<double(const double&)> transport_network_interpolation_function = transport_network_interpolation_function_parsing(argm["transport-network-interpolation-function"].back());
	int community_size = std::stoi(argm["community-size"].back());
	std::vector<int> community_degree = community_network_degree_parsing(argm["community-degree"].back());
	if (community_degree.size() > 1 && community_size != community_degree.size())
	{
		community_size = community_degree.size();
	}

	double mobility_rate = std::stod(argm["mobility-rate"].back());
	double community_infection_rate = std::stod(argm["community-infection-rate"].back());
	double transport_infection_rate = std::stod(argm["transport-infection-rate"].back());
	double recovery_rate = std::stod(argm["recovery-rate"].back());
	double immunity_loss_rate = std::stod(argm["immunity-loss-rate"].back());

	double initial_site = std::stoi(argm["initial-site"].back());
	double initial_prevalence = std::stod(argm["initial-prevalence"].back());

	double fractional_exponent = std::stod(argm["fractional-exponent"].back());

	double time = std::stod(argm["time"].back());

	std::string output_file = argm["output-file"].back();

	bool verbose = (bool)std::stoi(argm["verbose"].back());


	
	epidemic_transport_model _epidemic_transport_model;

	
	switch (transport_network_files.size())
	{
		case 1:
		{
			igraph_t transport_network;
			igraph_read_graph_graphmlfile(&transport_network, transport_network_files[0]);

			_epidemic_transport_model = epidemic_transport_model(transport_network, community_size, community_degree, initial_site, initial_prevalence, fractional_exponent, mobility_rate, community_infection_rate, transport_infection_rate, recovery_rate, immunity_loss_rate);

			igraph_cattribute_remove_all(&transport_network, true, true, true);
			igraph_destroy(&transport_network);
			break;
		}

		case 2:
		{
			std::array<igraph_t,2> transport_networks;
			for (size_t r = 0; r < transport_networks.size(); r++)
			{
				igraph_read_graph_graphmlfile(&transport_networks[r], transport_network_files[r]);
			}

			_epidemic_transport_model = epidemic_transport_model(transport_networks, transport_network_interpolation_function, community_size, community_degree, initial_site, initial_prevalence, fractional_exponent, mobility_rate, community_infection_rate, transport_infection_rate, recovery_rate, immunity_loss_rate);

			for (size_t r = 0; r < transport_networks.size(); r++)
			{
				igraph_cattribute_remove_all(&transport_networks[r], true, true, true);
				igraph_destroy(&transport_networks[r]);
			}
			break;
		}

		default:
		{
			break;
		}
	}


	_epidemic_transport_model.simulate(time);


	if (verbose)
	{
		std::cout << _epidemic_transport_model << std::endl;
	}

	if (strcmp(output_file.c_str(), "") != 0)
	{
		std::ofstream fstream;

		fstream.open(output_file);
		fstream << _epidemic_transport_model << std::endl;
		fstream.close();
	}
	else
	{
		if (!verbose)
		{
			std::cout << _epidemic_transport_model << std::endl;
		}
	}
	

	auto end_time = std::chrono::steady_clock::now();
	if (verbose)
	{
		std::chrono::duration<double> elapsed_time = end_time - start_time;
		std::cout << "*** END main (" << elapsed_time.count() << "s) ***" << std::endl;
	}

	return 0;
}
