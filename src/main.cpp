#include <iostream>
#include <fstream>

#include <getopt.h>
#include <map>

#include "epidemic_transport_model.hpp"

std::vector<double> linspace(double min, double max, size_t N)
{
	std::vector<double> range;

	double delta = (max - min) / double(N);
	for (int i = 0; i < (N + 1); i++) {
		range.push_back(min + i * delta);
	}
	return range;
}

std::map<std::string, std::string> argparse(int argc, char** argv)
{
	std::map<std::string, std::string> argm { {"transport-network-file", ""}, {"community-size", ""}, {"community-degree", ""}, {"prevalence", ""}, {"mobility-rate", ""}, {"community-infection-rate", ""}, {"transport-infection-rate", ""}, {"recovery-rate", ""}, {"immunity-loss-rate", ""}, {"time", ""}, {"output-file", ""}, {"verbose", "0"}, };

	const char* const short_opts = "w:N:k:p:m:B:b:g:s:T:o:vh";
	const option long_opts[] = {
		{"transport-network-file", required_argument, nullptr, 'w'},
		{"community-size", required_argument, nullptr, 'N'},
		{"community-degree", required_argument, nullptr, 'k'},
		{"prevalence", required_argument, nullptr, 'p'},
		{"mobility-rate", required_argument, nullptr, 'm'},
		{"community-infection-rate", required_argument, nullptr, 'B'},
		{"transport-infection-rate", required_argument, nullptr, 'b'},
		{"recovery-rate", required_argument, nullptr, 'g'},
		{"immunity-loss-rate", required_argument, nullptr, 's'},
		{"time", required_argument, nullptr, 'T'},
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
				argm["transport-network-file"] = optarg;
				break;
			case 'N':
				argm["community-size"] = optarg;
				break;
			case 'k':
				argm["community-degree"] = optarg;
				break;
			case 'p':
				argm["prevalence"] = optarg;
				break;
			case 'm':
				argm["mobility-rate"] = optarg;
				break;
			case 'B':
				argm["community-infection-rate"] = optarg;
				break;
			case 'b':
				argm["transport-infection-rate"] = optarg;
				break;
			case 'g':
				argm["recovery-rate"] = optarg;
				break;
			case 's':
				argm["immunity-loss-rate"] = optarg;
				break;
			case 'T':
				argm["time"] = optarg;
				break;
			case 'o':
				argm["output-file"] = optarg;
				break;
			case 'v':
				argm["verbose"] = "1";
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


	std::map<std::string, std::string> argm = argparse(argc, argv);

	std::string transport_network_file = argm["transport-network-file"];
	int community_size = std::stoi(argm["community-size"]);
	int community_degree = std::stoi(argm["community-degree"]);
	double prevalence = std::stod(argm["prevalence"]);
	double mobility_rate = std::stod(argm["mobility-rate"]);
	double community_infection_rate = std::stod(argm["community-infection-rate"]);
	double transport_infection_rate = std::stod(argm["transport-infection-rate"]);
	double recovery_rate = std::stod(argm["recovery-rate"]);
	double immunity_loss_rate = std::stod(argm["immunity-loss-rate"]);

	double time = std::stod(argm["time"]);

	std::string output_file = argm["output-file"];

	bool verbose = (bool)std::stoi(argm["verbose"]);

	igraph_t transport_network;

	FILE *ifile;

	ifile = fopen(transport_network_file.c_str(), "r");
	if (ifile == 0) {
		std::cout << "Cannot open transport network file." << std::endl;
	}

	igraph_read_graph_graphml(&transport_network, ifile, 0);
	fclose(ifile);


	epidemic_transport_model _epidemic_transport_model(&transport_network, community_size, community_degree, prevalence, mobility_rate, community_infection_rate, transport_infection_rate, recovery_rate, immunity_loss_rate);

	if (verbose)
	{
		std::cout << _epidemic_transport_model << std::endl;
	}

	_epidemic_transport_model.simulate(time);
	
	std::ofstream fstream;
	fstream.open(output_file);
	fstream << _epidemic_transport_model << std::endl;
	fstream.close();


	igraph_cattribute_remove_all(&transport_network, true, true, true);
	igraph_destroy(&transport_network);
	
	if (verbose)
	{
		std::cout << "*** END main ***" << std::endl;
	}

	return 0;
}
