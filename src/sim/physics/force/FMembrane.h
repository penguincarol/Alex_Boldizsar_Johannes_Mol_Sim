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

    void operator()() override;

    void setParticleContainer(ParticleContainer& pc) override;

    pair_fun_t& getForceFunction() override;

    fpair_fun_t getFastForceFunction() override;
};

} // sim::physics::force