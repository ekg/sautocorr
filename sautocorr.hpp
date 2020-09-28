#pragma once

#include <vector>
#include <numeric>
#include <cmath>
#include <iostream>
#include <cstdint>
#include <algorithm>

namespace sautocorr {

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
                       uint64_t k,
                       uint64_t stride,
                       double mean_i, double mean_j,
                       double stdev_i, double stdev_j) {
    //double mean_i = vec_mean(v.begin(), v.end()-k);
    //double stdev_i = vec_stdev(v.begin(), v.end()-k, mean_i);
    //double mean_j = vec_mean(v.begin()+k, v.end());
    //double stdev_j = vec_stdev(v.begin()+k, v.end(), mean_j);
    double sum_v = 0;
    for (size_t i = 0 ; i < v.size() - k ; i+=stride) {
        sum_v += (v[i] - mean_i) * (v[i+k] - mean_j);
    }
    return (sum_v / ((v.size() - k)/stride)) / (stdev_i * stdev_j);
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

struct repeat_t {
    uint64_t length = 0;
    double z_score = 0;
};

repeat_t repeat(const std::vector<uint8_t>& vals,
                uint64_t min_lag,
                uint64_t max_lag,
                uint64_t min_repeat,
                double target_z,
                uint64_t stride,
                const std::string& seq_name = "") {

    // bound max lag
    max_lag = std::min(max_lag, vals.size()/2);
    std::vector<double> autocorrs(max_lag - min_lag);
    double mean = vec_mean(vals.begin(), vals.end());
    double stdev = vec_stdev(vals.begin(), vals.end(), mean);
    for (int k = min_lag; k < max_lag; ++k) {
        autocorrs[k-min_lag] = autocorrelation(vals, k, stride, mean, mean, stdev, stdev);
    }

    double mean_ac = vec_mean(autocorrs.begin(), autocorrs.end());
    double stdev_ac = vec_stdev(autocorrs.begin(), autocorrs.end(), mean_ac);
    //std::cerr << "mean = " << mean_ac << std::endl;
    //std::cerr << "stdev = " << stdev_ac << std::endl;
    if (std::isnan(mean_ac) || std::isnan(stdev_ac)) {
        //std::cerr << "singularity found" << std::endl;
        return {0, 0};
    }

    std::vector<double> zscores(autocorrs.size());
    std::transform(autocorrs.begin(), autocorrs.end(),
                   zscores.begin(), [&](const double& d) {
                                        return (d - mean_ac)/stdev_ac;
                                    });

    if (seq_name.size()) {
        //std::cout << "seq.name\tlag\tautocorr\tz.score" << std::endl;
        for (int i = 0; i < autocorrs.size(); ++i) {
            std::cout << seq_name << "\t"
                      << i+min_lag << "\t"
                      << autocorrs[i] << "\t"
                      << zscores[i] << std::endl;
        }
    }

    // find our likely first max length
    double max_z = 0;
    uint64_t idx = 0;
    // find the peak
    // we expect it to be the highest score other than low
    bool seen_target = false;
    for (uint64_t i = 0; i < zscores.size(); ++i) {
        auto& z = zscores[i];
        if (z > target_z) {
            seen_target = true;
            max_z = std::max(z, max_z);
            idx = i;
        } else if (seen_target) {
            break;
        }
    }
    if (idx == 0) {
        //std::cerr << "no repeat found" << std::endl;
        return {0, 0};
    }
    uint64_t repeat_size = idx + min_lag;
    //std::cerr << "repeat is " << repeat_size << std::endl;
    return {repeat_size, max_z};
    // pick a repeat length
    // scan forward with it

    /*
    std::cout << "pos\tautocorr" << std::endl;
    // find the starts
    for (uint64_t i = 0; i < vals.size() - repeat_size * 2; ++i) {
        double a = autocorrelation_at(vals, repeat_size+10, i, repeat_size / 10);
        std::cout << i << "\t" << a << std::endl;
    }
    */

}

}
