//
// Created by alex on 11.01.23.
//

#include "FMembraneOMP.h"

namespace sim::physics::force {

    void FMembraneOMP::operator()() {
        // TODO impl me
    }

    void FMembraneOMP::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
    }

    pair_fun_t &FMembraneOMP::getForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};;
    }

    fpair_fun_t FMembraneOMP::getFastForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};
    }

    fpair_fun_alt_t FMembraneOMP::getFastForceAltFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};
    }
}
