//
// Created by alex on 11.01.23.
//

#include "FMembranePull.h"

namespace sim::physics::force {
    void FMembranePull::operator()() {
        particleContainer.runOnMembranes([&](std::vector<Membrane>& membranes,
                                            std::vector<double>& force,
                                            std::vector<double>& x,
                                            unsigned long count){
            for(auto& mem : membranes) {
                if (mem.getPullEndTime() < current_time) continue;

                for(auto& point : mem.getPullIndices()) {
                    size_t id = mem.getMembrNodes()[point[0]][point[1]];
                    force[id*3 + 0] += mem.getPullForce()[0];
                    force[id*3 + 1] += mem.getPullForce()[1];
                    force[id*3 + 2] += mem.getPullForce()[2];
                }
            }
        });
        current_time += delta_t;
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