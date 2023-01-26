//
// Created by alex on 26.11.2022.
//

#include "XStoermerVelvet.h"
#include <iostream>

namespace sim::physics::position {
    void XStoermerVelvet::operator()() {
#ifdef slow
        particleContainer.forAllParticles([&](Particle &p) {
            Eigen::Vector3d x = delta_t * p.getV() + delta_t * delta_t * p.getOldF() / (2 * p.getM());
            p.add_to_X(x);
        });
#else
        particleContainer.runOnActiveData([&](
                Kokkos::View<double*> &force,
                Kokkos::View<double*> &oldForce,
                Kokkos::View<double*> &x,
                Kokkos::View<double*> &v,
                Kokkos::View<double*> &m,
                Kokkos::View<int*> &type,
                unsigned long count,
                Kokkos::View<double*> &eps,
                Kokkos::View<double*> &sig,
                std::unordered_map<unsigned long, unsigned long> &id_to_index,
                std::vector<unsigned long> &activeParticles){
            for (unsigned long & activeParticle : activeParticles) {
                unsigned long index = id_to_index[activeParticle];
                x[3*index + 0] += delta_t * v[3*index + 0] + delta_t * delta_t * oldForce[3*index + 0] / (2 * m[index]);
                x[3*index + 1] += delta_t * v[3*index + 1] + delta_t * delta_t * oldForce[3*index + 1] / (2 * m[index]);
                x[3*index + 2] += delta_t * v[3*index + 2] + delta_t * delta_t * oldForce[3*index + 2] / (2 * m[index]);
            }
        });
#endif
    }
} //sim::physics::position
