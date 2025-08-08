#ifndef SCHRODINGER_SOLVER
#define SCHRODINGER_SOLVER

#include <Eigen/Eigen> // Requires Eigen library installation


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
        r_min = 1e-8; // Avoid singularity at r=0
        r_max = 50.0; // Bohr radii units
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

        // Solve eigenvalue problem
        Eigen::SparseMatrix<double> I(n_points, n_points);
        I.setIdentity();
        Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<double>> solver;
        solver.compute(H, Eigen::EigenvaluesOnly);
        
        // Extract ground state (lowest energy)
        int ground_state = 0;
        wavefunction = solver.eigenvectors().col(ground_state);
        normalize();
    }

    void normalize() {
        // Simpson's rule integration for radial probability
        double norm = 0.0;
        for (int i = 1; i < n_points; ++i) {
            double f1 = std::pow(wavefunction[i], 2);
            double f2 = std::pow(wavefunction[i-1], 2);
            norm += (f1 + f2) * dr;
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
};

#endif // SCHRODINGER_SOLVER