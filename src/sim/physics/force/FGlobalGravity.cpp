//
// Created by alex on 06.12.2022.
//

#include "FGlobalGravity.h"

void sim::physics::force::FGlobalGravity::operator()() {
    //perform gravity addition
    // we do not care if a particle is still active or not, faster this way
    particleContainer.runOnActiveData([&](std::vector<double> &force,
                                   std::vector<double> &oldForce,
                                   std::vector<double> &x,
                                   std::vector<double> &v,
                                   std::vector<double> &m,
                                   std::vector<int> &type,
                                   unsigned long count,
                                   std::vector<double> &eps,
                                   std::vector<double> &sig,
                                   std::unordered_map<unsigned long, unsigned long> &id_to_index,
                                   std::vector<unsigned long> &activeParticles) {
        for (auto [_,index] : id_to_index) {
            force[index*3 + 0] += m[index] * gGrav0;
            force[index*3 + 1] += m[index] * gGrav1;
            force[index*3 + 2] += m[index] * gGrav2;
        }
    });
}

void sim::physics::force::FGlobalGravity::setParticleContainer(ParticleContainer &pc) {
    particleContainer = pc;
}
