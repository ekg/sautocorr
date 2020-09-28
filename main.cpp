#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "seqiter.hpp"
#include "sautocorr.hpp"

int main (int argc, char** argv) {

    if (argc != 6) {
        std::cerr << argv[0] << "<sequence file> <min repeat> <max repeat> <min z-score> <stride>" << std::endl
                  << std::endl
                  << "    sequence autocorrelation detector" << std::endl
                  << std::endl
                  << "Input sequence file in FASTA or FASTQ format." << std::endl
                  << "max_lag scales the running time quadratically." << std::endl;
        return -1;
    }

    // load the input file
    std::string infile(argv[1]);
    int min_repeat = std::stoi(argv[2]);
    int max_repeat = std::stoi(argv[3]);
    double min_z = std::stod(argv[4]);
    int stride = std::stoi(argv[5]);

    std::cout << "seq.name\tlag\tautocorr\tz.score" << std::endl;
    seqiter::for_each_seq_in_file(
        infile,
        [&](const std::string& name,
            const std::string& seq) {
            std::vector<uint8_t> vec(seq.begin(), seq.end());
            sautocorr::repeat_t result = sautocorr::repeat(vec, name, min_repeat,
                                                           max_repeat, min_repeat, min_z, stride);
            std::cerr << name << "\t" << result.length << "\t" << result.z_score << std::endl;
        });

    return 0;
}
