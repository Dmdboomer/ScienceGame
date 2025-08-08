#include "ElemParticles.h"
#include "NonElemSubatomic.h"
#include "gnuplot-iostream.h"
#include "SchrodingerSolver.h"
#include "RadialDistribution.h"

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>


using namespace std;

class HydrogenAtom {
private:
    Proton* proton;
    Electron* electron;
    mt19937 gen;

public:
    HydrogenAtom() : proton(new Proton()), electron(new Electron()) {
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

    double waveFunction(double r) const {
        // ψ(r) = (1/√π) e⁻ʳ, where r is in units of Bohr radius
        return (1.0 / std::sqrt(M_PI)) * std::exp(-r);
    }

    double probabilityDensity(double r) const {
        double psi = waveFunction(r);
        return psi * psi; // |ψ(r)|²
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
            double prob = 4 * M_PI * r * r * probabilityDensity(r); // Radial probability
            data.push_back({r, prob});
        }
        return data;
    }

    // Generate 3D orbital slice data (xy-plane)
    void generate3DData(const std::string& filename, double size = 5.0, int resolution = 50) {
        std::ofstream outfile(filename);
        double step = 2 * size / resolution;
        
        for(double x = -size; x <= size; x += step) {
            for(double y = -size; y <= size; y += step) {
                double r = std::sqrt(x*x + y*y);
                double density = probabilityDensity(r);
                outfile << x << " " << y << " " << density << "\n";
            }
            outfile << "\n";
        }
    }

    //Actual 3d visual 
    void generateRandomPoints(const string& filename, int numPoints = 10000) {
        ofstream outfile(filename);
        uniform_real_distribution<> dis01(0.0, 1.0);
        
        for (int i = 0; i < numPoints; i++) {
            // Radial distance sampling (quantum mechanical distribution)
            double u1 = max(dis01(gen), 1e-15);
            double u2 = max(dis01(gen), 1e-15);
            double u3 = max(dis01(gen), 1e-15);
            double r = -0.5 * log(u1 * u2 * u3);
            
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

    // 2d
    /*
    Gnuplot gp;
    gp << "set title 'Hydrogen 1s Radial Probability'\n";
    gp << "set xlabel 'r (Bohr Radius)'\n";
    gp << "set ylabel 'Probability Density'\n";
    gp << "plot '-' with lines title '1s orbital'\n";
    gp.send1d(atom.radialDensityData());
    
    // Keep program open to view radial plot
    std::cin.get();
    return 0;
    */

    // Generate random points in the 1s orbital
    atom.generateRandomPoints("point_cloud.dat");
    
    // Plot points using Gnuplot with proper visualization settings
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