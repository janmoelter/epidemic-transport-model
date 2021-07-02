#include <iostream>
#include <fstream>

#include <vector>
#include <map>

#include <getopt.h>

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
	std::map<std::string, std::string> argm { {"transport-network-file", ""}, {"community-size", ""}, {"community-degree", ""}, {"prevalence", ""}, {"mobility-rate", ""}, {"community-infection-rate", ""}, {"transport-infection-rate", ""}, {"recovery-rate", ""}, {"time", ""}, {"output-file", ""}, };

	const char* const short_opts = "w:N:k:p:m:B:b:g:T:o:h";
	const option long_opts[] = {
		{"transport-network-file", required_argument, nullptr, 'w'},
		{"community-size", required_argument, nullptr, 'N'},
		{"community-degree", required_argument, nullptr, 'k'},
		{"prevalence", required_argument, nullptr, 'p'},
		{"mobility-rate", required_argument, nullptr, 'm'},
		{"community-infection-rate", required_argument, nullptr, 'B'},
		{"transport-infection-rate", required_argument, nullptr, 'b'},
		{"recovery-rate", required_argument, nullptr, 'g'},
		{"time", required_argument, nullptr, 'T'},
		{"output-file", required_argument, nullptr, 'o'},
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
			case 'T':
				argm["time"] = optarg;
				break;
			case 'o':
				argm["output-file"] = optarg;
				break;
			case 'h':
			case '?':
			default:
				break;
		}
	}

	return argm;
}

// ./bin/simulation -w transport_network.graphml -N 200 -k 5 -p 0.5 -m 10 -B 0.25 -b 0.25 -g 1 -T 10

int main(int argc, char **argv)
{
	std::map<std::string, std::string> argm = argparse(argc, argv);

	std::string transport_network_file = argm["transport-network-file"];
	int community_size = std::stoi(argm["community-size"]);
	int community_degree = std::stoi(argm["community-degree"]);
	double prevalence = std::stod(argm["prevalence"]);
	double mobility_rate = std::stod(argm["mobility-rate"]);
	double community_infection_rate = std::stod(argm["community-infection-rate"]);
	double transport_infection_rate = std::stod(argm["transport-infection-rate"]);
	double recovery_rate = std::stod(argm["recovery-rate"]);

	double time = std::stod(argm["time"]);

	std::string output_file = argm["output-file"];

	std::cout << std::setfill('#') << std::setw(80) << "" << std::endl;
	std::cout << "# transport network : " << transport_network_file << std::endl;
	std::cout << "# community size : " << community_size << std::endl;
	std::cout << "# community degree : " << community_degree << std::endl;
	std::cout << "# prevalence : " << prevalence << std::endl;
	std::cout << "# mobility rate : " << mobility_rate << std::endl;
	std::cout << "# community infection rate : " << community_infection_rate << std::endl;
	std::cout << "# transport infection rate : " << transport_infection_rate << std::endl;
	std::cout << "# recovery rate : " << recovery_rate << std::endl;
	std::cout << std::setfill('#') << std::setw(80) << "" << std::endl;
	std::cout << std::setfill(' ');

	epidemic_transport_model _epidemic_transport_model(transport_network_file, community_size, community_degree, prevalence, mobility_rate, community_infection_rate, transport_infection_rate, recovery_rate);
	std::list<std::tuple<double, int, double>> timeseries = _epidemic_transport_model.simulate(time);


	std::ofstream fstream;
	fstream.open(output_file);
	
	fstream << std::setfill('#') << std::setw(80) << "" << std::endl;
	fstream << "# transport network : " << transport_network_file << std::endl;
	fstream << "# community size : " << community_size << std::endl;
	fstream << "# community degree : " << community_degree << std::endl;
	fstream << "# prevalence : " << prevalence << std::endl;
	fstream << "# mobility rate : " << mobility_rate << std::endl;
	fstream << "# community infection rate : " << community_infection_rate << std::endl;
	fstream << "# transport infection rate : " << transport_infection_rate << std::endl;
	fstream << "# recovery rate : " << recovery_rate << std::endl;
	fstream << std::setfill('#') << std::setw(80) << "" << std::endl;
	fstream << std::setfill(' ');
	
	for (std::tuple<double, int, double> const& e: timeseries)
	{
		fstream << std::fixed << std::setprecision(9) << std::get<0>(e) << '\t' << std::get<1>(e) << '\t' << std::setprecision(12) << std::get<2>(e) << std::endl;
	}
	fstream.close();
	
	std::cout << "*** END main ***" << std::endl;

	return 0;
}
