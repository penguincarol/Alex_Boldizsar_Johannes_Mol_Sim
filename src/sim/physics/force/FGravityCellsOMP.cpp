//
// Created by alex on 11.01.23.
//

#include "FGravityCellsOMP.h"

namespace sim::physics::force {
    void FGravityCellsOMP::setPairFun() {
        pairFun = forceDelegate.getForceFunction();
        fpairFun = forceDelegate.getFastForceFunction();
    }

    void FGravityCellsOMP::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
        forceDelegate.setParticleContainer(pc);
        setPairFun();
    }

    pair_fun_t &FGravityCellsOMP::getForceFunction() {
        return pairFun;
    }

    fpair_fun_t FGravityCellsOMP::getFastForceFunction() {
        return fpairFun;
    }

    void FGravityCellsOMP::operator()() {
        // TODO write me
    }
} // force