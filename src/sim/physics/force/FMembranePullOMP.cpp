//
// Created by alex on 11.01.23.
//

#include "FMembranePullOMP.h"

namespace sim::physics::force {
    void FMembranePullOMP::operator()() {
        particleContainer.runOnMembranes([&](std::vector<Membrane>& membranes,
                                             std::vector<double>& force,
                                             std::vector<double>& x,
                                             unsigned long count,
                                             std::unordered_map<unsigned long, unsigned long> &id_to_index){
            for(auto& mem : membranes) {
                if (mem.getPullEndTime() < current_time) continue;

                #pragma omp parallel for default(none) shared(mem, force, id_to_index)
                for(auto& point : mem.getPullIndices()) {
                    size_t id = mem.getMembrNodes()[point[0]][point[1]];
                    size_t index = id_to_index[id];
                    force[index * 3 + 0] += mem.getPullForce()[0];
                    force[index * 3 + 1] += mem.getPullForce()[1];
                    force[index * 3 + 2] += mem.getPullForce()[2];
                }
            }
        });
        current_time += delta_t;
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