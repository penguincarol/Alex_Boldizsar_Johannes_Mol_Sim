//
// Created by alex on 26.11.2022.
//

#include "XStoermerVelvetOMP.h"

namespace sim::physics::position {
    void XStoermerVelvetOMP::operator()() {
        double delta_t = this->delta_t;
        particleContainer.runOnActiveData([delta_t](Kokkos::View<double*> &force,
                                       Kokkos::View<double*> &oldForce,
                                       Kokkos::View<double*> &x,
                                       Kokkos::View<double*> &v,
                                       Kokkos::View<double*> &m,
                                       auto,
                                       unsigned long count,
                                       auto, auto,std::unordered_map<unsigned long, unsigned long> &id_to_index,
                                       std::vector<unsigned long> &activeParticles) {
            unsigned long index;
#pragma omp parallel default(none) shared(force, oldForce, x, v, m, count, delta_t, activeParticles, id_to_index) private(index)
            {
#pragma omp for
                for (unsigned long & activeParticle : activeParticles) {
                    index = id_to_index[activeParticle];
                    x[index*3 + 0] += delta_t * v[index*3 + 0] + delta_t * delta_t * oldForce[index*3 + 0] / (2 * m[index]);
                    x[index*3 + 1] += delta_t * v[index*3 + 1] + delta_t * delta_t * oldForce[index*3 + 1] / (2 * m[index]);
                    x[index*3 + 2] += delta_t * v[index*3 + 2] + delta_t * delta_t * oldForce[index*3 + 2] / (2 * m[index]);
                }
            }
        });
    }
} //sim::physics::position