//
// Created by alex on 26.11.2022.
//

#include "VStoermerVelvetOMP.h"

namespace sim::physics::velocity {
    void VStoermerVelvetOMP::operator()() {
        double delta_t = this->delta_t;
        particleContainer.runOnData([&](std::vector<Particle>& particles,
                                        std::vector<Membrane>& membranes,
                                        ParticleContainer::VectorCoordWrapper& cells,
                                        unsigned long count,
                                        std::vector<unsigned long>& activeParticles,
                                        std::unordered_map<unsigned long, unsigned long> &id_to_index){
            unsigned long index;
            #pragma omp parallel default(none) shared(particles, delta_t, activeParticles, id_to_index) private(index)
            {
                #pragma omp for
                for (unsigned long i : activeParticles) {
                    index = id_to_index[i];
                    Particle& p = particles[index];
                    p.add_to_V(delta_t * (p.getOldF() + p.getF()) / (2 * p.getM()));
                }
            }
        });
    }
} //sim::physics::velocity