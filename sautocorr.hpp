#include <vector>
#include <numeric>
#include <cmath>

inline double vec_sum(const std::vector<double> v) {
    return std::accumulate(v.begin(), v.end(), 0.0);
}

inline double vec_mean(const std::vector<double>& v) {
    return vec_sum(v) / v.size();
}

inline double vec_stdev(const std::vector<double>& v) {
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double mean = vec_mean(v);
    return std::sqrt(sq_sum / v.size() - mean * mean);
}

template<typename T>
double autocorrelation(const std::vector<T>& v, uint64_t k) {
    double result = 0;
    for (size_t i = 0 ; i < v.size() - k ; i++) {
        result += std::abs(v[i] - v[i+k]);
    }
    return (double)result / (double)(v.size() - k);
}

template<typename T>
double autocorrelation_at(const std::vector<T>& v, uint64_t offset, uint64_t k) {
    double result = 0;
    for (size_t i = offset ; i < offset + k ; i++) {
        result += std::abs(v[i] - v[i+k]);
    }
    return (double)result / (double)k;
}
