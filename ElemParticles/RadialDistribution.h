#ifndef RADIAL_DISTRIBUTION
#define RADIAL_DISTRIBUTION

#include <Eigen/Eigen>
#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <random>
#include <algorithm>

class RadialDistribution {
private:
    std::vector<double> r_values;
    std::vector<double> cdf;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;

public:
    RadialDistribution(const Eigen::VectorXd& r_vec, 
                        const Eigen::VectorXd& psi,
                        int points) {
        r_values.resize(points);
        cdf.resize(points);
        
        // Build CDF: ∫|ψ(r)|² r² dr (radial probability)
        double integral = 0.0;
        for (int i = 0; i < psi.size(); ++i) {
            double r = r_vec[i];
            double prob = std::pow(psi[i], 2) * r * r;
            if (i > 0) {
                integral += 0.5 * (prob + prev_prob) * dr;
            }
            cdf[i] = integral;
            r_values[i] = r;
            prev_prob = prob;
        }
        // Normalize CDF to [0,1]
        cdf /= integral; 

        // RNG setup
        std::random_device rd;
        gen.seed(rd());
        dist = std::uniform_real_distribution(0.0, 1.0);
    }

    double sample() {
        double u = dist(gen);
        auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
        int idx = std::distance(cdf.begin(), it);
        return r_values[idx]; // Inverse transform sampling
    }
};

#endif // RADIAL_DISTRIBUTION