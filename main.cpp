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

    uint64_t max_j = 5000;

    for (uint64_t i = 0; i < vals.size(); ++i) {
        for (uint64_t j = 10; j < max_j; j+=500) {
            range_t r = next_repeat_range(vals, i, j, max_lag);
            std::cerr << i << "\t" << j << "\t" << r.begin << "\t" << r.end << std::endl;
            if (r.end > 0) {
                i = r.end;
                break;
            }
        }
    }
                              

    return 0;
}
