//
// Created by alex on 11.01.23.
//

#include "FMembranePullOMP.h"

namespace sim::physics::force {
    void FMembranePullOMP::operator()() {
        //TODO impl
    }

    void FMembranePullOMP::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
    }

    pair_fun_t &FMembranePullOMP::getForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};;
    }

    fpair_fun_t FMembranePullOMP::getFastForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};;
    }
} // force