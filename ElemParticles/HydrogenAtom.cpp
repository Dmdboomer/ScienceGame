#include "ElemParticles.h"
#include "NonElemSubatomic.h"
#include "gnuplot-iostream.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <Eigen/Eigen>

using namespace std;

//  SchrodingerSolver with public getters
class SchrodingerSolver {
private:
    double reduced_mass;
    double alpha = 7.2973525693e-3; // Fine structure constant
    int n_points;
    double r_min, r_max, dr;
    Eigen::VectorXd r_vector;
    Eigen::VectorXd wavefunction;

public:
    SchrodingerSolver(double mu, int points = 1000)
        : reduced_mass(mu), n_points(points) {
        double bohr_radius = 1.0 / (reduced_mass * alpha);
        r_min = 1e-8 * bohr_radius;
        r_max = 10.0 * bohr_radius; // Cover up to 10a₀
        dr = (r_max - r_min) / (n_points - 1);
        solve();
}

    void solve() {
        // Create radial grid
        r_vector.resize(n_points);
        for (int i = 0; i < n_points; ++i) {
            r_vector[i] = r_min + i * dr;
        }

        // Build Hamiltonian matrix (tridiagonal)
        Eigen::SparseMatrix<double> H(n_points, n_points);
        std::vector<Eigen::Triplet<double>> triplets;
        
        double hbar2 = 1.0; // ħ=1 in natural units
        double ke_prefactor = -hbar2 / (2.0 * reduced_mass * dr * dr);
        double diag_factor = -2.0 * ke_prefactor;

        for (int i = 1; i < n_points - 1; ++i) {
            // Potential term: V(r) = -α/r
            double V = -alpha / r_vector[i];
            
            // Central diagonal
            triplets.emplace_back(i, i, diag_factor + V);
            
            // Adjacent diagonals
            triplets.emplace_back(i, i-1, ke_prefactor);
            triplets.emplace_back(i, i+1, ke_prefactor);
        }
        H.setFromTriplets(triplets.begin(), triplets.end());

        // Solve eigenvalue problem (WITH eigenvectors)
        Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<double>> solver;
        solver.compute(H, Eigen::ComputeEigenvectors); // Key fix
        if (solver.info() != Eigen::Success) {
            std::cerr << "Eigensolver failed." << std::endl;
            return;
        }

        // Extract ground state (now eigenvectors are accessible)
        wavefunction = solver.eigenvectors().col(0);
        normalize();
    }

    void normalize() {
        double norm = 0.0;
        // Corrected trapezoidal integration
        for (int i = 1; i < n_points; ++i) {
            double f_prev = std::pow(wavefunction[i-1], 2);
            double f_current = std::pow(wavefunction[i], 2);
            norm += (f_prev + f_current) * dr / 2.0; // Factor of 1/2 added
        }
        wavefunction /= std::sqrt(norm);
    }

    double radialWavefunction(double r) const {
        if (r < r_min || r > r_max) return 0.0;
        // Linear interpolation
        double index = (r - r_min) / dr;
        int idx = static_cast<int>(index);
        double frac = index - idx;
        return wavefunction[idx] * (1.0 - frac) + wavefunction[idx+1] * frac;
    }

    const Eigen::VectorXd& getRVector() const { return r_vector; }
    const Eigen::VectorXd& getWavefunction() const { return wavefunction; }
    int getNumPoints() const { return n_points; }

};

// Fixed RadialDistribution with correct integral calculation
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
        // Assuming uniform grid
        double dr = r_vec[1] - r_vec[0];
        r_values.resize(points);
        cdf.resize(points);
        double integral = 0.0;
        double prev_prob = 0.0;

        for (int i = 0; i < points; ++i) {
            double r = r_vec[i];
            // |u(r)|² is the radial probability density
            double prob = std::pow(psi[i], 2);

            if (i == 0) {
                integral = 0.0;
            } else {
                // Trapezoidal rule integration
                integral += 0.5 * (prev_prob + prob) * dr;
            }

            r_values[i] = r;
            cdf[i] = integral;
            prev_prob = prob;
        }

        // Normalize CDF to [0,1]
        double normalization = cdf[points-1];
        for (int i = 0; i < points; ++i) {
            cdf[i] /= normalization;
        }

        // Initialize random number generator
        std::random_device rd;
        gen.seed(rd());
        dist = std::uniform_real_distribution<double>(0.0, 1.0);
    }

    double sample() {
        double u = dist(gen);
        // Find the smallest index such that cdf[i] >= u
        auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
        int idx = std::distance(cdf.begin(), it);
        // Clamp index to valid range
        idx = std::min(idx, static_cast<int>(r_values.size()) - 1);
        return r_values[idx];
    }
};

class HydrogenAtom {
private:
    Proton* proton;
    Electron* electron;
    SchrodingerSolver solver;
    RadialDistribution radialDist;
    mt19937 gen;

public:
    HydrogenAtom() : 
        proton(new Proton()), 
        electron(new Electron()),
        solver(reducedMass()),
        radialDist(solver.getRVector(), solver.getWavefunction(), solver.getNumPoints())
    {
        random_device rd;
        gen.seed(rd());
    }

    ~HydrogenAtom() {
        delete proton;
        delete electron;
    }

    double reducedMass() const {
        double m_e = electron->getMass();
        double m_p = proton->getMass();
        return (m_e * m_p) / (m_e + m_p);
    }

    double groundStateEnergy() const {
        const double ALPHA = 7.2973525693e-3; // Fine structure constant
        double mu = reducedMass(); // Reduced mass in GeV/c²
        return -0.5 * mu * ALPHA * ALPHA; // E = -0.5 μ α² (GeV)
    }

    // Radial wavefunction u(r) = r*R(r)
    double waveFunction(double r) const {
        return solver.radialWavefunction(r);
    }

    // Radial probability density |u(r)|²
    double probabilityDensity(double r) const {
        double psi = waveFunction(r);
        return psi * psi;
    }

    double getBohrRadius() const {
        const double ALPHA = 7.2973525693e-3;
        double mu = reducedMass();
        return 1.0 / (mu * ALPHA); // a₀ = 1/(μ α) (GeV⁻¹)
    }

    // Generate radial probability density data
    std::vector<std::pair<double, double>> radialDensityData(double max_r = 10.0, int points = 100) {
        std::vector<std::pair<double, double>> data;
        for(int i = 0; i <= points; ++i) {
            double r = max_r * i / points;
            // Use |u(r)|² directly
            double prob = probabilityDensity(r);
            data.push_back({r, prob});
        }
        return data;
    }

    //Actual 3d visual 
    void generateRandomPoints(const string& filename, int numPoints = 10000) {
        ofstream outfile(filename);
        uniform_real_distribution<> dis01(0.0, 1.0);
        
        for (int i = 0; i < numPoints; i++) {
            // Radial distance sampling (using numerical radial distribution)
            double r = radialDist.sample();
            
            // Angular distribution (uniform sphere)
            double phi = 2 * M_PI * dis01(gen);
            double theta = acos(2.0 * dis01(gen) - 1.0);
            
            // Convert to Cartesian coordinates
            double x = r * sin(theta) * cos(phi);
            double y = r * sin(theta) * sin(phi);
            double z = r * cos(theta);
            
            outfile << x << " " << y << " " << z << "\n";
        }
    }

};

int main() {
    HydrogenAtom atom;

    // Generate random points in the 1s orbital
    atom.generateRandomPoints("point_cloud.dat");
    
    // Plot points using Gnuplot
    Gnuplot gp;
    gp << "set title 'Hydrogen 1s Orbital Electron Probability Distribution'\n";
    gp << "set xlabel 'x (Bohr Radii)'\n";
    gp << "set ylabel 'y (Bohr Radii)'\n";
    gp << "set zlabel 'z (Bohr Radii)'\n";
    gp << "set view equal xyz\n";
    gp << "set style line 1 lc rgb '#0066FF' pt 7 ps 0.1\n";
    gp << "splot 'point_cloud.dat' with points ls 1 title 'Electron Position'\n";
    
    // Keep program open
    cout << "Plot displayed. Press Enter to exit.\n";
    cin.get();
    return 0;
}