//
// Created by alex on 11.01.23.
//

#include "FGlobalGravityOMP.h"
#include "defaults.h"

namespace sim::physics::force {

    void FGlobalGravityOMP::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
    }

    pair_fun_t &FGlobalGravityOMP::getForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};
    }

    fpair_fun_t FGlobalGravityOMP::getFastForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};
    }

    void FGlobalGravityOMP::operator()() {
        // TODO impl
    }
} // force