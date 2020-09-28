#include "sautocorr.hpp"

namespace sautocorr {

repeat_t repeat(const std::vector<uint8_t>& vals,
                uint64_t min_lag,
                uint64_t max_lag,
                uint64_t min_repeat,
                double target_z,
                uint64_t stride,
                const std::string& seq_name) {

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

