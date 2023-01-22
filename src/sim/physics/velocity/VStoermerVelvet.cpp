//
// Created by alex on 26.11.2022.
//

#include "VStoermerVelvet.h"

namespace sim::physics::velocity {
    void VStoermerVelvet::operator()() {
#ifdef slow
        particleContainer.forAllParticles([&](Particle &p) {
            Eigen::Vector3d v = delta_t * (p.getOldF() + p.getF()) / (2 * p.getM());
            p.add_to_V(v);
        });
#else
        particleContainer.runOnActiveData([&](
                vec4d_t &force,
                vec4d_t &oldForce,
                vec4d_t &x,
                vec4d_t &v,
                std::vector<double> &m, std::vector<int> &type,
                unsigned long count, std::vector<double> &eps, std::vector<double> &sig,
                std::unordered_map<unsigned long, unsigned long> &id_to_index,
                std::vector<unsigned long> &activeParticles){
            for (auto [_,index] : id_to_index) {
                v[4*index + 0] += delta_t * (oldForce[4*index + 0] + force[4*index + 0]) / (2 * m[index]);
                v[4*index + 1] += delta_t * (oldForce[4*index + 1] + force[4*index + 1]) / (2 * m[index]);
                v[4*index + 2] += delta_t * (oldForce[4*index + 2] + force[4*index + 2]) / (2 * m[index]);
            }
        });
#endif
    }
} //sim::physics::velocity