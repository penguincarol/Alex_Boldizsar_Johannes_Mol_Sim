//
// Created by alex on 11.01.23.
//

#include "FGlobalGravityOMP.h"
#include "defaults.h"

namespace sim::physics::force {

    void FGlobalGravityOMP::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
    }

    void FGlobalGravityOMP::operator()() {
        particleContainer.runOnActiveData([&](vec4d_t &force,
                                              vec4d_t &oldForce,
                                              vec4d_t &x,
                                              vec4d_t &v,
                                              std::vector<double> &m,
                                              std::vector<int> &type,
                                              unsigned long count,
                                              std::vector<double> &eps,
                                              std::vector<double> &sig,
                                              std::unordered_map<unsigned long, unsigned long> &id_to_index,
                                              std::vector<unsigned long> &activeParticles) {
            //generate tasks
            std::vector<std::vector<unsigned long>> tasks;
            tasks.emplace_back();
            for (auto[_,index] : id_to_index) {
                if(tasks.back().size() > max_thread_tasks) tasks.emplace_back();
                tasks.back().emplace_back(index);
            }

            double gGrav0 = this->gGrav0;
            double gGrav1 = this->gGrav1;
            double gGrav2 = this->gGrav2;
            #pragma omp parallel for default(none) shared(tasks, force, m, gGrav0, gGrav1, gGrav2)
            for (auto& task : tasks) {
                for (unsigned long index : task) {
                    force[index*4 + 0] += m[index] * gGrav0;
                    force[index*4 + 1] += m[index] * gGrav1;
                    force[index*4 + 2] += m[index] * gGrav2;
                }
            }
        });
    }
} // force