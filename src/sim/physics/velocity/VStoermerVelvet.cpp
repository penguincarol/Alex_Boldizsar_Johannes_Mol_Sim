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
        particleContainer.runOnActiveData([&](Kokkos::View<double*> &force, Kokkos::View<double*> &oldForce,
                                        Kokkos::View<double*> &x, Kokkos::View<double*> &v,
                                        Kokkos::View<double*> &m, Kokkos::View<int*> &type,
                                        unsigned long count, Kokkos::View<double*> &eps, Kokkos::View<double*> &sig,
                                        std::unordered_map<unsigned long, unsigned long> &id_to_index,
                                        std::vector<unsigned long> &activeParticles){
            for (auto [_,index] : id_to_index) {
                v[3*index + 0] += delta_t * (oldForce[3*index + 0] + force[3*index + 0]) / (2 * m[index]);
                v[3*index + 1] += delta_t * (oldForce[3*index + 1] + force[3*index + 1]) / (2 * m[index]);
                v[3*index + 2] += delta_t * (oldForce[3*index + 2] + force[3*index + 2]) / (2 * m[index]);
            }
        });
#endif
    }
} //sim::physics::velocity