#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "sautocorr.hpp"

int main (int argc, char** argv) {

    if (argc != 3) {
        std::cerr << argv[0] << "<sequence file> <max lag>" << std::endl
                  << std::endl
                  << "    sequence autocorrelation detector" << std::endl
                  << std::endl
                  << "Input sequence file should contain a single integer per line." << std::endl
                  << "The max lag to calculate bounds running time quadratically." << std::endl;
        return -1;
    }

    // load the input file
    std::string infile(argv[1]);
    int max_lag = std::stoi(argv[2]);

    std::ifstream f(infile.c_str());
    std::vector<int64_t> vals;
    if (f.is_open()) {
        std::string line;
        while (f >> line) {
            vals.push_back(std::atoi(line.c_str()));
        }
    } else {
        std::cerr << "Unable to open file:  " << infile << std::endl;
        return -1;
    }

    std::vector<double> autocorrs;
    for (int k = 1; k < std::min((int)max_lag, (int)vals.size()); ++k) {
        autocorrs.push_back(autocorrelation(vals, k));
    }

    std::cout << "lag\tautocorr" << std::endl;
    for (int i = 0; i < autocorrs.size(); ++i) {
        std::cout << i+1 << "\t" << autocorrs[i] << std::endl;
    }

    return 0;
}
