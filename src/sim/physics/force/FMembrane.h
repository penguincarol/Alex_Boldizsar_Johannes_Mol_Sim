//
// Created by johnny on 09.01.23.
//
#pragma once

#include "ForceFunctorBase.h"

namespace sim::physics::force {
/**
 * calculate the force caused by springs in every membrane
 * */
class FMembrane : public ForceFunctorBase {
private:
    pair_fun_t pairFun;
    fpair_fun_t fpairFun;
    ForceFunctorBase* forceDelegate;
    std::vector<Membrane> membranes;

    void setPairFun();

public:
    /**
     * the created instance will take ownership of ff and will delete it upon deconstruction.
     * */
    FMembrane(
                       double st,
                       double et,
                       double dt,
                       double eps,
                       double sig,
                       ParticleContainer &pc,
                       ForceFunctorBase* ff
    ) : ForceFunctorBase(st, et, dt, eps, sig, pc), forceDelegate(ff), membranes(pc.getMembranes()) {
        setPairFun();
        forceDelegate->setParticleContainer(pc);
    }

    ~FMembrane() override {
        delete forceDelegate;
    }

    /**
     * Helper method to increase readability of the code
     * Not meant to be used on its own
     * computes the force between the particles p1 and p2 that are in the same membrane caused by the spring between them
     */
    void addSpringForce(size_t p1i, size_t p1j, size_t p2i, size_t p2j,
                                           Membrane& membrane, std::vector<double>& force, std::vector<double>& x);

    void operator()() override;

    void setParticleContainer(ParticleContainer& pc) override;

    pair_fun_t& getForceFunction() override;

    fpair_fun_t getFastForceFunction() override;
};

} // sim::physics::force