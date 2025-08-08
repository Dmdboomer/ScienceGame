#ifndef NON_ELEM_SUBATOMIC_H
#define NON_ELEM_SUBATOMIC_H

#include "ElemParticles.h"

class Proton : public Fermion {
public:
    Proton(bool anti = false)
        : Fermion(0.938, 1.0, anti, "Proton") {
        adjustForAntiparticle();
    }

    ElementaryParticle* createAntiparticle() const override {
        return new Proton(true);
    }

    virtual ~Proton() = default;
};


class Neutron : public Fermion {
public:
    Neutron(bool anti = false)
        : Fermion(0.939, 0.0, anti, "Neutron") {
        adjustForAntiparticle();
    }

    ElementaryParticle* createAntiparticle() const override {
        return new Neutron(true);
    }

    virtual ~Neutron() = default;
};

#endif // NON_ELEM_SUBATOMIC_H