#include <iostream>
#include <fstream>

#include <getopt.h>
#include <map>

#include <regex>

#include "epidemic_transport_model.hpp"


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

std::vector<double> linspace(double min, double max, size_t N)
{
	std::vector<double> range;

	double delta = (max - min) / double(N);
	for (int i = 0; i < (N + 1); i++) {
		range.push_back(min + i * delta);
	}
	return range;
}

std::map<std::string, std::vector<std::string>> argparse(int argc, char** argv)
{
	std::map<std::string, std::vector<std::string>> argm { {"transport-network-file", {}}, {"transport-network-interpolation-function", {}}, {"community-size", {}}, {"community-degree", {}}, {"mobility-rate", {}}, {"community-infection-rate", {}}, {"transport-infection-rate", {}}, {"recovery-rate", {}}, {"immunity-loss-rate", {}}, {"time", {}}, {"initial-site", {"inf"}}, {"initial-prevalence", {}}, {"output-file", {}}, {"verbose", {"0"}}, };

	const char* const short_opts = "w:f:N:k:m:B:b:g:s:T:x:p:o:vh";
	const option long_opts[] = {
		{"transport-network-file", required_argument, nullptr, 'w'},
		{"transport-network-interpolation-function", required_argument, nullptr, 'w'},
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
		{"output-file", required_argument, nullptr, 'o'},
		{"verbose", no_argument, nullptr, 'v'},
		{"help", no_argument, nullptr, 'h'},
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
			case 'o':
				argm["output-file"].push_back(optarg);
				break;
			case 'v':
				argm["verbose"].push_back("1");
				break;
			case 'h':
			case '?':
			default:
				break;
		}
	}

	return argm;
}

int main(int argc, char **argv)
{
	igraph_set_attribute_table(&igraph_cattribute_table);

	std::cout << "+ " << __LINE__ << std::endl;

	std::map<std::string, std::vector<std::string>> argm = argparse(argc, argv);

	std::vector<std::string> transport_network_files = argm["transport-network-file"];
	std::string transport_network_interpolation_function_string = argm["transport-network-interpolation-function"].back();
	int community_size = std::stoi(argm["community-size"].back());
	int community_degree = std::stoi(argm["community-degree"].back());
	double mobility_rate = std::stod(argm["mobility-rate"].back());
	double community_infection_rate = std::stod(argm["community-infection-rate"].back());
	double transport_infection_rate = std::stod(argm["transport-infection-rate"].back());
	double recovery_rate = std::stod(argm["recovery-rate"].back());
	double immunity_loss_rate = std::stod(argm["immunity-loss-rate"].back());

	double initial_site = std::stoi(argm["initial-site"].back());
	double initial_prevalence = std::stod(argm["initial-prevalence"].back());

	double time = std::stod(argm["time"].back());

	std::string output_file = argm["output-file"].back();

	bool verbose = (bool)std::stoi(argm["verbose"].back());

	std::cout << "+ " << __LINE__ << std::endl;

	std::array<igraph_t,2> transport_networks;
	switch (transport_network_files.size())
	{
		case 1:
		{
			igraph_read_graph_graphmlfile(&transport_networks[0], transport_network_files[0]);
			igraph_copy(&transport_networks[1], &transport_networks[0]);
			break;
		}

		case 2:
		{
			igraph_read_graph_graphmlfile(&transport_networks[0], transport_network_files[0]);
			igraph_read_graph_graphmlfile(&transport_networks[1], transport_network_files[1]);
			break;
		}

		default:
		{
			break;
		}
	}

	std::array<igraph_t *,2> _transport_networks = {&transport_networks[0], &transport_networks[1]};


	std::function<double(const double&)> transport_network_interpolation_function;

	std::regex function_string_regex("^(\\w*)\\(([+]?\\d*(?:\\.(?:[0-9]+))?);([+-]?\\d*(?:\\.(?:[0-9]+))?)\\)\\+\\+([+]?\\d*(?:\\.(?:[0-9]+))?)$");
	std::smatch matches;
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

		transport_network_interpolation_function = [f, arg_offset, f_width, f_period] (double t)->double { return f(fmin(fmod(t - (arg_offset - f_width/2.), f_period) / f_width, 1.)) * (t > arg_offset - f_width/2.); };
	}
	else
	{
		throw std::invalid_argument(transport_network_interpolation_function_string);
	}


	epidemic_transport_model _epidemic_transport_model(_transport_networks, transport_network_interpolation_function, community_size, community_degree, initial_site, initial_prevalence, mobility_rate, community_infection_rate, transport_infection_rate, recovery_rate, immunity_loss_rate);

	
	if (verbose)
	{
		std::cout << _epidemic_transport_model << std::endl;
	}
	
	_epidemic_transport_model.simulate(time);
	
	std::ofstream fstream;
	fstream.open(output_file);
	fstream << _epidemic_transport_model << std::endl;
	fstream.close();


	igraph_cattribute_remove_all(&transport_networks[0], true, true, true);
	igraph_cattribute_remove_all(&transport_networks[1], true, true, true);

	igraph_destroy(&transport_networks[0]);
	igraph_destroy(&transport_networks[1]);
	
	if (verbose)
	{
		std::cout << "*** END main ***" << std::endl;
	}
	

	return 0;
}
