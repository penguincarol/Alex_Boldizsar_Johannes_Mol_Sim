//
// Created by alex on 26.11.2022.
//

#include "VStoermerVelvetOMP.h"

namespace sim::physics::velocity {
    void VStoermerVelvetOMP::operator()() {
        double delta_t = this->delta_t;
        particleContainer.runOnActiveData([delta_t](vec4d_t &force,
                                                    vec4d_t &oldForce,
                                                    vec4d_t &x,
                                                    vec4d_t &v,
                                                    std::vector<double> &m,
                                                    std::vector<int> &type,
                                                    unsigned long count, auto, auto,
                                                    std::unordered_map<unsigned long, unsigned long> &id_to_index,
                                                    std::vector<unsigned long> &activeParticles) {
            unsigned long index;
#pragma omp parallel default(none) shared(force, oldForce, v, m, count, delta_t, activeParticles, id_to_index) private(index)
            {
#pragma omp for
                for (unsigned long i : activeParticles) {
                    index = id_to_index[i];
                    v[index*4 + 0] += delta_t * (oldForce[index*4 + 0] + force[index*4 + 0]) / (2 * m[index]);
                    v[index*4 + 1] += delta_t * (oldForce[index*4 + 1] + force[index*4 + 1]) / (2 * m[index]);
                    v[index*4 + 2] += delta_t * (oldForce[index*4 + 2] + force[index*4 + 2]) / (2 * m[index]);
                }
            }
        });
    }
} //sim::physics::velocity