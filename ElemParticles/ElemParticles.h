#ifndef ELEMENTARY_PARTICLES_H
#define ELEMENTARY_PARTICLES_H

#include <string>
#include <stdexcept>
#include <cmath>

// Base class for all elementary particles
class ElementaryParticle {
protected:
    double mass;        // Mass in GeV/c²
    double charge;      // Charge in elementary charge units
    double spin;        // Spin in ħ units
    bool isAntiparticle;
    std::string name;

    // Set properties for antiparticles
    void adjustForAntiparticle() {
        if (!isSelfConjugate() && isAntiparticle) {
            charge = -charge;
        }
    }

    // Check if particle is its own antiparticle
    bool isSelfConjugate() const {
        return (name == "Photon" || name == "Gluon" || 
                name == "ZBoson" || name == "HiggsBoson");
    }

public:
    ElementaryParticle(double m, double c, double s, bool anti, std::string n)
        : mass(m), charge(c), spin(s), isAntiparticle(anti), name(n) {}

    // Getters for particle properties
    double getMass() const { return mass; }
    double getCharge() const { return charge; }
    double getSpin() const { return spin; }
    std::string getName() const { 
        return (isAntiparticle && !isSelfConjugate()) ? "anti-" + name : name;
    }
    bool getIsAntiparticle() const { return isAntiparticle; }

    // Clone method for creating antiparticles
    virtual ElementaryParticle* createAntiparticle() const = 0;
    virtual ~ElementaryParticle() = default;
};

// Fermion base class (half-integer spin)
class Fermion : public ElementaryParticle {
public:
    Fermion(double m, double c, bool anti, std::string n)
        : ElementaryParticle(m, c, 0.5, anti, n) {}
};

// Boson base class (integer spin)
class Boson : public ElementaryParticle {
public:
    Boson(double m, double c, double s, bool anti, std::string n)
        : ElementaryParticle(m, c, s, anti, n) {}
};

// Lepton base class
class Lepton : public Fermion {
public:
    Lepton(double m, double c, bool anti, std::string n)
        : Fermion(m, c, anti, n) {}
};

// Quark base class
class Quark : public Fermion {
public:
    Quark(double m, double c, bool anti, std::string n)
        : Fermion(m, c, anti, n) {}
};

// Specific particle implementations
class Electron : public Lepton {
public:
    Electron(bool anti = false)
        : Lepton(0.000511, -1.0, anti, "Electron") {}

    ElementaryParticle* createAntiparticle() const override {
        return new Electron(true);
    }
};

class ElectronNeutrino : public Lepton {
public:
    ElectronNeutrino(bool anti = false)
        : Lepton(0.0, 0.0, anti, "ElectronNeutrino") {}

    ElementaryParticle* createAntiparticle() const override {
        return new ElectronNeutrino(true);
    }
};

class Muon : public Lepton {
public:
    Muon(bool anti = false)
        : Lepton(0.1057, -1.0, anti, "Muon") {}

    ElementaryParticle* createAntiparticle() const override {
        return new Muon(true);
    }
};

class MuonNeutrino : public Lepton {
public:
    MuonNeutrino(bool anti = false)
        : Lepton(0.0, 0.0, anti, "MuonNeutrino") {}

    ElementaryParticle* createAntiparticle() const override {
        return new MuonNeutrino(true);
    }
};

class Tau : public Lepton {
public:
    Tau(bool anti = false)
        : Lepton(1.777, -1.0, anti, "Tau") {}

    ElementaryParticle* createAntiparticle() const override {
        return new Tau(true);
    }
};

class TauNeutrino : public Lepton {
public:
    TauNeutrino(bool anti = false)
        : Lepton(0.018, 0.0, anti, "TauNeutrino") {}

    ElementaryParticle* createAntiparticle() const override {
        return new TauNeutrino(true);
    }
};

class UpQuark : public Quark {
public:
    UpQuark(bool anti = false)
        : Quark(0.0022, 2.0/3, anti, "UpQuark") {}

    ElementaryParticle* createAntiparticle() const override {
        return new UpQuark(true);
    }
};

class DownQuark : public Quark {
public:
    DownQuark(bool anti = false)
        : Quark(0.0047, -1.0/3, anti, "DownQuark") {}

    ElementaryParticle* createAntiparticle() const override {
        return new DownQuark(true);
    }
};

class CharmQuark : public Quark {
public:
    CharmQuark(bool anti = false)
        : Quark(1.28, 2.0/3, anti, "CharmQuark") {}

    ElementaryParticle* createAntiparticle() const override {
        return new CharmQuark(true);
    }
};

class StrangeQuark : public Quark {
public:
    StrangeQuark(bool anti = false)
        : Quark(0.096, -1.0/3, anti, "StrangeQuark") {}

    ElementaryParticle* createAntiparticle() const override {
        return new StrangeQuark(true);
    }
};

class TopQuark : public Quark {
public:
    TopQuark(bool anti = false)
        : Quark(173.1, 2.0/3, anti, "TopQuark") {}

    ElementaryParticle* createAntiparticle() const override {
        return new TopQuark(true);
    }
};

class BottomQuark : public Quark {
public:
    BottomQuark(bool anti = false)
        : Quark(4.18, -1.0/3, anti, "BottomQuark") {}

    ElementaryParticle* createAntiparticle() const override {
        return new BottomQuark(true);
    }
};

class Photon : public Boson {
public:
    Photon()
        : Boson(0.0, 0.0, 1.0, false, "Photon") {}

    ElementaryParticle* createAntiparticle() const override {
        return new Photon();  // Self-conjugate
    }
};

class Gluon : public Boson {
public:
    Gluon()
        : Boson(0.0, 0.0, 1.0, false, "Gluon") {}

    ElementaryParticle* createAntiparticle() const override {
        return new Gluon();  // Self-conjugate
    }
};

class WBoson : public Boson {
public:
    WBoson(bool isPositive = true)
        : Boson(80.4, isPositive ? 1.0 : -1.0, 1.0, !isPositive, "WBoson") {}

    ElementaryParticle* createAntiparticle() const override {
        return new WBoson(charge < 0);  // Flips charge state
    }
};

class ZBoson : public Boson {
public:
    ZBoson()
        : Boson(91.2, 0.0, 1.0, false, "ZBoson") {}

    ElementaryParticle* createAntiparticle() const override {
        return new ZBoson();  // Self-conjugate
    }
};

class HiggsBoson : public Boson {
public:
    HiggsBoson()
        : Boson(125.3, 0.0, 0.0, false, "HiggsBoson") {}

    ElementaryParticle* createAntiparticle() const override {
        return new HiggsBoson();  // Self-conjugate
    }
};

#endif // ELEMENTARY_PARTICLES_H