#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "autocorr.hpp"

int main (int argc, char** argv) {

    // load the input file
    std::string infile(argv[1]);
    int max_lag = std::stoi(argv[2]);
    std::ifstream f(infile.c_str());
    std::stringstream buffer;
    buffer << f.rdbuf();
    std::string str = buffer.str(); // hmm copy
    //std::cerr << str;

    std::cout << "lag\tautocorr" << std::endl;
    for (int k = 1; k < std::min((int)max_lag, (int)str.size()); ++k) {
        std::cout << k << "\t" << autocorrelation(str, k) << std::endl;
    }

    return 0;
}
