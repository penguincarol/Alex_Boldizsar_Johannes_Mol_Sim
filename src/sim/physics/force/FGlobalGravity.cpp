//
// Created by alex on 06.12.2022.
//

#include "FGlobalGravity.h"

void sim::physics::force::FGlobalGravity::operator()() {
    //perform gravity addition
    // we do not care if a particle is still active or not, faster this way
    particleContainer.runOnActiveData([&](Kokkos::View<double*> &force,
                                   Kokkos::View<double*> &oldForce,
                                   Kokkos::View<double*> &x,
                                   Kokkos::View<double*> &v,
                                   Kokkos::View<double*> &m,
                                   Kokkos::View<int*> &type,
                                   unsigned long count,
                                   Kokkos::View<double*> &eps,
                                   Kokkos::View<double*> &sig,
                                   std::unordered_map<unsigned long, unsigned long> &id_to_index,
                                   std::vector<unsigned long> &activeParticles) {
        for (auto [_,index] : id_to_index) {
            force[index*3 + 0] += std::max(m[index], 0.0) * gGrav0;
            force[index*3 + 1] += std::max(m[index], 0.0) * gGrav1;
            force[index*3 + 2] += std::max(m[index], 0.0) * gGrav2;
        }
    });
}

void sim::physics::force::FGlobalGravity::setParticleContainer(ParticleContainer &pc) {
    particleContainer = pc;
}
