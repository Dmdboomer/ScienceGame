#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

const double GRAVITATIONAL_CONSTANT = 6.67430e-11;

struct GravityFormula {
    int type; // 0 = exponent, 1 = distance offset
    string expression;
    double modifier;
};

GravityFormula generate_formula() {
    int formula_type = rand() % 1;
    GravityFormula formula;
    formula.type = formula_type;

    formula.modifier = (rand() % 11) / 10.0;
    stringstream modifier_stream;
    modifier_stream << fixed << setprecision(1) << formula.modifier;
    string mod_str = modifier_stream.str();

    switch(formula_type) {
        case 0:
            formula.expression = "G * m1 * m2 / pow(r + " + mod_str + ", 2)";
            break;
        case 1:
            double exponent = 2.0 + formula.modifier;
            stringstream exponent_stream;
            exponent_stream << fixed << setprecision(1) << exponent;
            formula.expression = "G * m1 * m2 / pow(r, " + exponent_stream.str() + ")";
            break;
    }
    
    return formula;
}

double calculate_force(GravityFormula formula, double m1, double m2, double r) {
    switch(formula.type) {
        case 0: 
            return GRAVITATIONAL_CONSTANT * m1 * m2 / pow(r, 2 + formula.modifier);
        case 1:
            return GRAVITATIONAL_CONSTANT * m1 * m2 / pow(r + formula.modifier, 2);
        default:
            return 0;
    }
}

int main() {
    srand(time(0));

    // Newton's apple experiment
    GravityFormula custom_formula = generate_formula();
    double m_earth = 5.972e24;
    double m_apple = 0.2;
    cout << "Newton's Apple Experiment (Modified Formula):\n";
    cout << "\n The mass of the apple: " 
        <<  setprecision(3) << m_apple << "\n";
    cout << "--------------------------------------------------\n";
    cout << "Height (m) | Calculated Acceleration (m/s²)\n";
    cout << "--------------------------------------------------\n";
    
    for(double height = 1.0; height <= 10.0; height += 1.0) {
        double r = 1e4 + height;
        double force = calculate_force(custom_formula, m_earth, m_apple, r);
        double acceleration = force / m_apple;
        
        cout << fixed << setprecision(1) << setw(9) << height << " | "
             << scientific << setprecision(4) << acceleration << endl;
    }
    
    // User challenge
    string user_input;
    cout << "\nEnter the ACTUAL formula used:\n> ";
    cin.ignore();
    getline(cin, user_input);
    
    if(user_input.find(custom_formula.expression.substr(0, 25)) != string::npos) {
        cout << "\nCORRECT! Actual formula: " << custom_formula.expression << "\n\n";
    } else {
        cout << "\nINCORRECT! Actual formula: " << custom_formula.expression << "\n\n";
    }
    
    // Torsion balance experiment (FIXED: Added mass/distance details)
    double G_TORSION = 5.0 + (rand() % 1001) / 20.0; // 5.00 - 55.00
    double mass1 = 0.01, mass2 = 0.01; // 10g each
    double distance = 0.05; // 5cm
    
    cout << "Cavendish Torsion Balance Experiment:\n";
    cout << "--------------------------------------------------\n";
    cout << "Small Mass 1: " << mass1 << " kg\n";
    cout << "Small Mass 2: " << mass2 << " kg\n";
    cout << "Distance: " << distance << " m\n";
    cout << "Measured Force: " << scientific << GRAVITATIONAL_CONSTANT * mass1 * mass2 / pow(distance, 2) << " N\n\n";
    
    // User challenge
    double user_G;
    cout << "Enter the gravitational constant used (XX.XXe-XX format): ";
    cin >> user_G;
    
    if(fabs(user_G - G_TORSION) < 0.005) {
        cout << "\nCORRECT! G = " << fixed << setprecision(2) << G_TORSION << endl;
    } else {
        cout << "\nINCORRECT! Actual G = " << fixed << setprecision(2) << G_TORSION << endl;
    }
    
    return 0;
}