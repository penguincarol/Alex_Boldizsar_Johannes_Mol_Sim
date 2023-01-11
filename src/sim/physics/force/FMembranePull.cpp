//
// Created by alex on 11.01.23.
//

#include "FMembranePull.h"

namespace sim::physics::force {
    void FMembranePull::operator()() {
        //TODO impl me
    }

    void FMembranePull::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
    }

    pair_fun_t &FMembranePull::getForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};;;
    }

    fpair_fun_t FMembranePull::getFastForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};;
    }
} // force