//
// Created by alex on 11.01.23.
//

#include "FMembrane.h"

namespace sim::physics::force {
    void FMembrane::operator()() {
        // TODO impl me
    }

    void FMembrane::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
    }

    pair_fun_t &FMembrane::getForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};;
    }

    fpair_fun_t FMembrane::getFastForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};
    }
} // force