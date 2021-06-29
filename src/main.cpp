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
	std::map<std::string, std::string> argm { {"transport-network-file", ""}, {"community-size", ""}, {"community-degree", ""}, {"prevalence", ""}, {"mobility-rate", ""}, {"infection-rate", ""}, {"recovery-rate", ""}, {"coupling", ""}, {"time", ""}, {"output-file", ""}, };

	const char* const short_opts = "w:N:k:p:m:b:g:c:T:o:h";
	const option long_opts[] = {
		{"transport-network-file", required_argument, nullptr, 'w'},
		{"community-size", required_argument, nullptr, 'N'},
		{"community-degree", required_argument, nullptr, 'k'},
		{"prevalence", required_argument, nullptr, 'p'},
		{"mobility-rate", required_argument, nullptr, 'm'},
		{"infection-rate", required_argument, nullptr, 'b'},
		{"recovery-rate", required_argument, nullptr, 'g'},
		{"coupling", required_argument, nullptr, 'c'},
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
			case 'b':
				argm["infection-rate"] = optarg;
				break;
			case 'g':
				argm["recovery-rate"] = optarg;
				break;
			case 'c':
				argm["coupling"] = optarg;
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

// ./bin/simulation -w transport_network.graphml -N 200 -k 5 -p 0.5 -m 10 -b 0.25 -g 1 -c 1 -T 10

int main(int argc, char **argv)
{
	std::map<std::string, std::string> argm = argparse(argc, argv);

	std::string transport_network_file = argm["transport-network-file"];
	int community_size = std::stoi(argm["community-size"]);
	int community_degree = std::stoi(argm["community-degree"]);
	double prevalence = std::stod(argm["prevalence"]);
	double mobility_rate = std::stod(argm["mobility-rate"]);
	double infection_rate = std::stod(argm["infection-rate"]);
	double recovery_rate = std::stod(argm["recovery-rate"]);
	double coupling = std::stod(argm["coupling"]);

	double time = std::stod(argm["time"]);

	std::string output_file = argm["output-file"];

	std::cout << std::setfill('#') << std::setw(80) << "" << std::endl;
	std::cout << "# transport network : " << transport_network_file << std::endl;
	std::cout << "# community size : " << community_size << std::endl;
	std::cout << "# community degree : " << community_degree << std::endl;
	std::cout << "# prevalence : " << prevalence << std::endl;
	std::cout << "# mobility rate : " << mobility_rate << std::endl;
	std::cout << "# infection rate : " << infection_rate << std::endl;
	std::cout << "# recovery rate : " << recovery_rate << std::endl;
	std::cout << "# coupling parameter : " << coupling << std::endl;
	std::cout << std::setfill('#') << std::setw(80) << "" << std::endl;
	std::cout << std::setfill(' ');

	epidemic_transport_model _epidemic_transport_model(transport_network_file, community_size, community_degree, prevalence, mobility_rate, infection_rate, recovery_rate, coupling);
	std::list<std::tuple<double, int, double>> timeseries = _epidemic_transport_model.run(time);


	std::ofstream fstream;
	fstream.open(output_file);
	
	fstream << std::setfill('#') << std::setw(80) << "" << std::endl;
	fstream << "# transport network : " << transport_network_file << std::endl;
	fstream << "# community size : " << community_size << std::endl;
	fstream << "# community degree : " << community_degree << std::endl;
	fstream << "# prevalence : " << prevalence << std::endl;
	fstream << "# mobility rate : " << mobility_rate << std::endl;
	fstream << "# infection rate : " << infection_rate << std::endl;
	fstream << "# recovery rate : " << recovery_rate << std::endl;
	fstream << "# coupling parameter : " << coupling << std::endl;
	fstream << std::setfill('#') << std::setw(80) << "" << std::endl;
	fstream << std::setfill(' ');
	
	for (std::tuple<double, int, double> const& e: timeseries)
	{
		fstream << std::fixed << std::setprecision(9) << std::get<0>(e) << '\t' << std::get<1>(e) << '\t' << std::setprecision(12) << std::get<2>(e) << std::endl;
	}
	fstream.close();
	
	std::cout << "*** END main ***" << std::endl;


	//int trials = 10;
	//
	//double _mobility_rate = 1;
	//double _infection_rate;
	//double _recovery_rate = 10;
	//double _coupling;
	//
	//std::vector<double> beta_gamma = linspace(0.05, 0.25, 50);
	//std::vector<double> eta = linspace(0, 1, 50);
	//
	//std::vector<double> prevalence;
	//
	//std::list<std::tuple<double, int, double>> timeseries;
	//
	//std::cout << std::setfill('#') << std::setw(80) << "" << std::endl;
	//std::cout << "# community size : " << _community_size << std::endl;
	//std::cout << "# k : " << k << std::endl;
	//std::cout << "# mobility rate : " << _mobility_rate << std::endl;
	//std::cout << "# infection rate : " << "*" << std::endl;
	//std::cout << "# recovery rate : " << _recovery_rate << std::endl;
	//std::cout << "# coupling parameter : " << "*" << std::endl;
	//std::cout << std::setfill('#') << std::setw(80) << "" << std::endl;
	//
	//std::cout << " β/γ \\  η  \t";
	//for (int j = 0; j < eta.size(); j++)
	//{
	//	std::cout << std::fixed << std::setprecision(9) << eta[j] << '\t' << std::flush;
	//}
	//std::cout << std::endl;
	//
	//for (int i = 0; i < beta_gamma.size(); i++)
	//{
	//	_infection_rate = _recovery_rate * beta_gamma[i];
	//
	//	std::cout << std::fixed << std::setprecision(9) << beta_gamma[i] << '\t' << std::flush;
	//
	//	for (int j = 0; j < eta.size(); j++)
	//	{
	//		_coupling = eta[j];
	//
	//		for (int t = 0; t < trials; t++)
	//		{
	//			epidemic_transport_model _epidemic_transport_model(_community_size, k, 0.1, _mobility_rate, _infection_rate, _recovery_rate, _coupling);
	//
	//			timeseries = _epidemic_transport_model.run(5);
	//			
	//			prevalence.push_back(std::get<2>(timeseries.back()));
	//			//std::cout << std::get<2>(timeseries.back()) << std::endl;
	//		}
	//
	//		std::cout << std::fixed << std::setprecision(9) << std::accumulate(prevalence.begin(), prevalence.end(), 0.0) / prevalence.size() << '\t' << std::flush;
	//
	//		prevalence.clear();
	//	}
	//	std::cout << std::endl;
	//
	//}

	return 0;
}
