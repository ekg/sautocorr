#include "sautocorr.hpp"

range_t next_repeat_range(const std::vector<int64_t>& vals,
                          uint64_t start,
                          uint64_t length,
                          uint64_t max_lag,
                          uint64_t min_repeat,
                          double target_z) {

    std::vector<double> autocorrs;
    for (int k = start; k < start+std::min((int)max_lag, (int)length); ++k) {
        autocorrs.push_back(autocorrelation(vals, k));
    }

    double mean_ac = vec_mean(autocorrs.begin(), autocorrs.end());
    double stdev_ac = vec_stdev(autocorrs.begin(), autocorrs.end(), mean_ac);
    std::cerr << "mean = " << mean_ac << std::endl;
    std::cerr << "stdev = " << stdev_ac << std::endl;
    if (std::isnan(mean_ac) || std::isnan(stdev_ac)) {
        std::cerr << "singularity found" << std::endl;
        return {0, 0};
    }

    /*
    std::cout << "lag\tautocorr\tz.score" << std::endl;
    for (int i = 0; i < autocorrs.size(); ++i) {
        std::cout << i+1 << "\t"
                  << autocorrs[i] << "\t"
                  << (autocorrs[i] - mean_ac)/stdev_ac << std::endl;
    }
    */
    std::vector<double> zscores(autocorrs.size());
    std::transform(autocorrs.begin(), autocorrs.end(),
                   zscores.begin(), [&](const double& d) {
                                        return (d - mean_ac)/stdev_ac;
                                    });

    // find our likely first max length
    double max_z = 0;
    uint64_t idx = 0;
    // find the peak
    // we expect it to be the highest score other than low
    for (uint64_t i = min_repeat; i < zscores.size(); ++i) {
        if (zscores[i] > target_z) {
            max_z = zscores[i];
            idx = i;
            break;
        }
    }
    if (idx == 0) {
        std::cerr << "no repeat found" << std::endl;
        return {0, 0};
    }
    uint64_t repeat_size = idx;
    std::cerr << "repeat is " << repeat_size << std::endl;
    return {start, repeat_size};
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
