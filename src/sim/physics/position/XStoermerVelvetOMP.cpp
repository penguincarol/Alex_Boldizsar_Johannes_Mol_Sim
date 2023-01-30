//
// Created by alex on 26.11.2022.
//

#include "XStoermerVelvetOMP.h"

namespace sim::physics::position {
    void XStoermerVelvetOMP::operator()() {
        double delta_t = this->delta_t;
        particleContainer.runOnActiveData([delta_t](std::vector<double> &force,
                                       std::vector<double> &oldForce,
                                       std::vector<double> &x,
                                       std::vector<double> &v,
                                       std::vector<double> &m,
                                       auto,
                                       unsigned long count,
                                       auto, auto,
                                       std::vector<unsigned long> &activeParticles) {
            unsigned long index;
#pragma omp parallel default(none) shared(force, oldForce, x, v, m, count, delta_t, activeParticles) private(index)
            {
#pragma omp for
                for (unsigned long & activeParticle : activeParticles) {
                    index = activeParticle;
                    x[index*3 + 0] += delta_t * v[index*3 + 0] + delta_t * delta_t * oldForce[index*3 + 0] / (2 * m[index]);
                    x[index*3 + 1] += delta_t * v[index*3 + 1] + delta_t * delta_t * oldForce[index*3 + 1] / (2 * m[index]);
                    x[index*3 + 2] += delta_t * v[index*3 + 2] + delta_t * delta_t * oldForce[index*3 + 2] / (2 * m[index]);
                }
            }
        });
    }
} //sim::physics::position