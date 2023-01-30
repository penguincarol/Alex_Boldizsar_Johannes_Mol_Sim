//
// Created by alex on 06.12.2022.
//

#pragma once

#include "BoundsFunctorBase.h"
#include "sim/physics/force/ForceHandler.h"

#include <omp.h>

namespace sim::physics::bounds {
    template<sim::physics::bounds::side S>
    class BoundsPeriodic : public BoundsFunctorBase<S> {
    private:
        force::ForceHandler& forceHandler;
        bool mirrorMinor;
        bool mirrorMajor;

    public:
        ~BoundsPeriodic() override = default;

        /**
         * @brief Creates a periodic bound with the specified parameters
         * @param st
         * @param et
         * @param dt
         * @param eps
         * @param sig
         * @param pc
         * @param ff
         * @param mMinor
         * @param mMajor
         */
        BoundsPeriodic(double st, double et, double dt, double eps, double sig, ParticleContainer &pc, force::ForceHandler& fh, bool mMinor, bool mMajor, bool eOMP)
                : BoundsFunctorBase<S>(st, et, dt, eps, sig, pc, eOMP), forceHandler(fh), mirrorMinor(mMinor), mirrorMajor(mMajor) {}

        /**Reflects particle upon nearing the border.
         * Halo construction and force calculation has to be called independently by calling handle Halo.
         * */
        void operator()() final {
            // bounds handler will clear halo cell entries prior to any operator() call of a periodic bound
            // move particles that went beyond the bounds to the other side
            std::unordered_set<unsigned long> indices;
            this->particleContainer.template popExternalParticles<S>(indices);
            this->particleContainer.template moveExternalParticles<S>(indices);
            indices.clear();
        }

        bool isPeriodic() final {return true;}

        void generateHalo() final {
            // create halo particles for particles in this side and store in PC (this is tho only place that has knowledge of particles at edges and corners)
            this->particleContainer.template storeBorderParticlesToHalo<S>(mirrorMinor, mirrorMajor);
        }

        /**Calculates the force for this halo side in periodic bounds.*/
        void calcHaloForce() final {
            if(!this->enableOMP) {
                this->particleContainer.template forAllPairsHaloSide<S>(forceHandler.getFastForceFunction());
                return;
            }

            this->particleContainer.template runOnDataCell([&](std::vector<double> &force,
                                                               std::vector<double> &oldForce,
                                                               std::vector<double> &x,
                                                               std::vector<double> &v,
                                                               std::vector<double> &m,
                                                               std::vector<int> &t,
                                                               unsigned long count,
                                                               ParticleContainer::VectorCoordWrapper &cells,
                                                               std::vector<double> &eps,
                                                               std::vector<double> &sig){
                auto tasks = this->particleContainer.template generateHaloBorderTasks<S>();
                auto retFun = this->forceHandler.getFastForceRetFunction();
                double* f = force.data();
                size_t size = force.size();

                #pragma omp parallel for default(none) shared(tasks, x, t, eps,m,force, sig, size, retFun) reduction(+:f[:size])
                for(size_t tid = 0; tid < static_cast<size_t>(omp_get_max_threads()); tid++) {
                    for(size_t hcId = 0; hcId < tasks[tid].size(); hcId++) {
                        auto hCell = std::get<0>(tasks[tid][hcId][0]);
                        for (size_t pairIndex = 0; pairIndex < tasks[tid][hcId].size(); pairIndex++) {
                            auto &task = tasks[tid][hcId][pairIndex];
                            //auto hCell = std::get<0>(task);
                            auto bCell = std::get<1>(task);
                            auto off0 = std::get<2>(task);
                            auto off1 = std::get<3>(task);
                            auto off2 = std::get<4>(task);

                            for (size_t indexH = 0; indexH < hCell->size(); indexH++) {
                                unsigned long indexJ = (*hCell)[indexH];
                                double x0 = x[3 * indexJ + 0] + off0;
                                double x1 = x[3 * indexJ + 1] + off1;
                                double x2 = x[3 * indexJ + 2] + off2;

                                for (size_t indexB = 0; indexB < bCell->size(); indexB++) {
                                    unsigned long indexI = (*bCell)[indexB];
                                    auto new_force = retFun(force,x,eps,sig,m,t,indexI,x0,x1,x2,eps[indexJ],sig[indexJ],0,t[indexJ]);

                                    f[indexI * 3 + 0] -= new_force[0];
                                    f[indexI * 3 + 1] -= new_force[1];
                                    f[indexI * 3 + 2] -= new_force[2];
                                    f[indexJ * 3 + 0] += new_force[0];
                                    f[indexJ * 3 + 1] += new_force[1];
                                    f[indexJ * 3 + 2] += new_force[2];
                                }
                            }
                        }
                        hCell->clear();
                    }
                }
            });
        }
    };
} // sim::physics::bounds