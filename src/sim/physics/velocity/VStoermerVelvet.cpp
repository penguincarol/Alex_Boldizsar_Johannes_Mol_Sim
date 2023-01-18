//
// Created by alex on 26.11.2022.
//

#include "VStoermerVelvet.h"

namespace sim::physics::velocity {
    void VStoermerVelvet::operator()() {
        particleContainer.forAllParticles([&](Particle &p) {
            p.add_to_V(delta_t * (p.getOldF() + p.getF()) / (2 * p.getM()));
        });
    }
} //sim::physics::velocity