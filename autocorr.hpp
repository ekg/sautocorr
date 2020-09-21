#include <vector>

//typedef std::vector<char> Sequence;

double autocorrelation(const std::string& seq, int k) {
    //cassert(k > 0 && k < (int) seq.size());
    double result = 0;
    for (size_t i = 0 ; i < seq.size() - k ; i++) {
        result += std::abs(seq[i] - seq[i+k]);
    }
    return result / (seq.size() - k);
}

