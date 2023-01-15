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
        particleContainer.runOnActiveData([&](std::vector<double> &force, std::vector<double> &oldForce,
                                       std::vector<double> &x, std::vector<double> &v,
                                       std::vector<double> &m, std::vector<int> &type,
                                       unsigned long count, std::vector<double> &eps, std::vector<double> &sig,
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
