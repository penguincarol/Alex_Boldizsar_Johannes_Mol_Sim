//
// Created by alex on 26.11.2022.
//

#include "XStoermerVelvetOMP.h"

namespace sim::physics::position {
    void XStoermerVelvetOMP::operator()() {
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
                for (unsigned long & activeParticle : activeParticles) {
                    index = id_to_index[activeParticle];
                    Particle& p = particles[index];
                    p.add_to_F(delta_t * p.getV() + delta_t*delta_t*p.getOldF()/(2*p.getM()));
                }
            }
        });
    }
} //sim::physics::position