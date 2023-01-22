//
// Created by alex on 26.11.2022.
//

#include "XStoermerVelvetOMP.h"

namespace sim::physics::position {
    void XStoermerVelvetOMP::operator()() {
        double delta_t = this->delta_t;
        particleContainer.runOnActiveData([delta_t](
                vec4d_t &force,
                vec4d_t &oldForce,
                vec4d_t &x,
                vec4d_t &v,
                std::vector<double> &m,
                auto,
                unsigned long count,
                auto, auto,std::unordered_map<unsigned long, unsigned long> &id_to_index,
                std::vector<unsigned long> &activeParticles
                ) {
            unsigned long index;
#pragma omp parallel default(none) shared(force, oldForce, x, v, m, count, delta_t, activeParticles, id_to_index) private(index)
            {
#pragma omp for
                for (unsigned long & activeParticle : activeParticles) {
                    index = id_to_index[activeParticle];
                    x[index*4 + 0] += delta_t * v[index*4 + 0] + delta_t * delta_t * oldForce[index*4 + 0] / (2 * m[index]);
                    x[index*4 + 1] += delta_t * v[index*4 + 1] + delta_t * delta_t * oldForce[index*4 + 1] / (2 * m[index]);
                    x[index*4 + 2] += delta_t * v[index*4 + 2] + delta_t * delta_t * oldForce[index*4 + 2] / (2 * m[index]);
                }
            }
        });
    }
} //sim::physics::position