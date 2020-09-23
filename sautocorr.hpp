#pragma once

#include <vector>
#include <numeric>
#include <cmath>
#include <iostream>
#include <cstdint>
#include <algorithm>

template<typename T>
double vec_sum(T v_begin,
               T v_end) {
    return std::accumulate(v_begin, v_end, 0.0);
}

template<typename T>
double vec_mean(T v_begin,
                T v_end) {
    return vec_sum(v_begin, v_end) / (v_end - v_begin);
}

template<typename T>
double vec_stdev(T v_begin,
                 T v_end,
                 const double& mean) {
    double sq_sum = std::inner_product(v_begin, v_end, v_begin, 0.0);
    return std::sqrt(sq_sum / (v_end - v_begin) - mean * mean);
}

template<typename T>
double autocorrelation(const std::vector<T>& v,
                       uint64_t k) {
    double mean_i = vec_mean(v.begin(), v.end()-k);
    double stdev_i = vec_stdev(v.begin(), v.end()-k, mean_i);
    double mean_j = vec_mean(v.begin()+k, v.end());
    double stdev_j = vec_stdev(v.begin()+k, v.end(), mean_j);
    double sum_v = 0;
    for (size_t i = 0 ; i < v.size() - k ; i++) {
        sum_v += (v[i] - mean_i) * (v[i+k] - mean_j);
    }
    return (sum_v / (v.size() - k)) / (stdev_i * stdev_j);
}

template<typename T>
double autocorrelation_at(const std::vector<T>& v, uint64_t k, uint64_t j, uint64_t m) {
    double mean_i = vec_mean(v.begin()+j, v.begin()+j+m);
    double stdev_i = vec_stdev(v.begin()+j, v.begin()+j+m, mean_i);
    double mean_j = vec_mean(v.begin()+j+k, v.begin()+(j+m+k));
    double stdev_j = vec_stdev(v.begin()+j+k, v.begin()+(j+m+k), mean_j);
    double sum_v = 0;
    for (size_t i = j ; i < j + m ; i++) {
        sum_v += (v[i] - mean_i) * (v[i+k] - mean_j);
    }
    return (sum_v / (v.size() - k)) / (stdev_i * stdev_j);

    //return autocorrelation_range(v, k, 0, v.size());
}

struct range_t {
    uint64_t begin = 0;
    uint64_t end = 0;
};


range_t next_repeat_range(const std::vector<int64_t>& vals,
                          uint64_t start,
                          uint64_t length,
                          uint64_t max_lag,
                          uint64_t min_repeat = 10,
                          double target_z = 3);
