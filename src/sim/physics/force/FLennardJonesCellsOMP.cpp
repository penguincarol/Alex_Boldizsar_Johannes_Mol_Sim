//
// Created by alex on 11.01.23.
//

#include "FLennardJonesCellsOMP.h"
#include "defaults.h"
#include "Kokkos_ScatterView.hpp"

#include <iostream>
#include <omp.h>

namespace sim::physics::force {
    /**
     * @brief Returns the force function used
     *
     * @return pair_fun_t&
     */
    pair_fun_t &FLennardJonesCellsOMP::getForceFunction() {
        return pairFun;
    }

    /**
     * @brief The name says it all
     *
     * @param pc
     */
    void FLennardJonesCellsOMP::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
        forceDelegate.setParticleContainer(pc);
        setPairFun();
    }

    void FLennardJonesCellsOMP::setPairFun() {
        pairFun = forceDelegate.getForceFunction();
        fpairFun = forceDelegate.getFastForceFunction();
        fpairFunAlt = forceDelegate.getFastForceAltFunction();
        fpairFunRet = forceDelegate.getFastForceRetFunction();
    }

    fpair_fun_ret_t FLennardJonesCellsOMP::getFastForceRetFunction() {
        return fpairFunRet;
    }

    fpair_fun_t FLennardJonesCellsOMP::getFastForceFunction() {
        return fpairFun;
    }

    fpair_fun_alt_t FLennardJonesCellsOMP::getFastForceAltFunction() {
        return fpairFunAlt;
    }
    static const double rt3_2 = std::pow(2, 1.0 / 3.0);

#ifndef ONE_DIMENSIONAL_TASKS
#ifndef THREE_DIMENSIONAL_TASKS
#ifndef TASK_ORIENTED_2D
    void FLennardJonesCellsOMP::operator()() {
        particleContainer.runOnDataCell([&](Kokkos::View<double*> &force,
                                            Kokkos::View<double*> &oldForce,
                                            Kokkos::View<double*> &x,
                                            Kokkos::View<double*> &v,
                                            Kokkos::View<double*> &m,
                                            Kokkos::View<int*> &t,
                                            unsigned long count,
                                            ParticleContainer::VectorCoordWrapper& cells,
                                            Kokkos::View<double*> &eps,
                                            Kokkos::View<double*> &sig) {
            std::vector<std::pair<int,int>> pairs;
            for(auto& cell: cells) {
                for(size_t i = 0; i < cell.size(); i++){
                    for(size_t j = i+1; j < cell.size(); j++){
                        pairs.emplace_back(cell[i], cell[j]);
                    }
                }
            }

            Kokkos::View<int*> indIs = Kokkos::View<int*>("indI", pairs.size());
            Kokkos::View<int*> indJs = Kokkos::View<int*>("indI", pairs.size());
            for(int i = 0; i < pairs.size(); i++) {
                indIs[i] = pairs[i].first;
                indJs[i] = pairs[i].second;
            }

            Kokkos::Experimental::ScatterView<double*> _f(force);
            Kokkos::parallel_for("FLJOMPSingle", pairs.size(), KOKKOS_LAMBDA (const int &i) {
                int indexI = indIs[i];
                int indexJ = indJs[i];
                double sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1;
                sigma = (sig[indexI] + sig[indexJ]) / 2;
                sigma2 = sigma * sigma;
                sigma6 = sigma2 * sigma2 * sigma2;
                epsilon = std::sqrt(eps[indexI] * eps[indexJ]); // TODO this can be cached
                d0 = x[indexI*3 + 0] - x[indexJ*3 + 0];
                d1 = x[indexI*3 + 1] - x[indexJ*3 + 1];
                d2 = x[indexI*3 + 2] - x[indexJ*3 + 2];
                dsqr = d0*d0 + d1*d1 + d2*d2;
                //check if is membrane -> need to skip attractive forces
                if (t[indexI] & 0x80000000 || t[indexJ] & 0x80000000) {
                    if (dsqr >= rt3_2 * sigma2) return;
                }

                l2NInvSquare = 1 / (dsqr);
                fac0 = 24 * epsilon * l2NInvSquare;
                l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;
                fac1_sum1 = sigma6 * l2NInvPow6;
                fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);

                auto access = _f.access();
                access[indexI*3 + 0] -= fac0 * fac1 * d0;
                access[indexI*3 + 1] -= fac0 * fac1 * d1;
                access[indexI*3 + 2] -= fac0 * fac1 * d2;
                access[indexJ*3 + 0] += fac0 * fac1 * d0;
                access[indexJ*3 + 1] += fac0 * fac1 * d1;
                access[indexJ*3 + 2] += fac0 * fac1 * d2;
            });
            Kokkos::fence();

            const std::vector<std::pair<unsigned long, unsigned long>>& alternativeTaskGroups = particleContainer.generateFlatModel();
            size_t size = alternativeTaskGroups.size();
            Kokkos::View<int*> indIs2 = Kokkos::View<int*>("indI", size);
            Kokkos::View<int*> indJs2 = Kokkos::View<int*>("indI", size);
            for(int i = 0; i < size; i++) {
                indIs2[i] = alternativeTaskGroups[i].first;
                indJs2[i] = alternativeTaskGroups[i].second;
            }

            Kokkos::Experimental::ScatterView<double*> _ff(force);
            Kokkos::parallel_for("FLJOMPPair", size, KOKKOS_LAMBDA (const int &i) {
                int indexI = indIs[i];
                int indexJ = indJs[i];
                double sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1;
                sigma = (sig[indexI] + sig[indexJ]) / 2;
                sigma2 = sigma * sigma;
                sigma6 = sigma2 * sigma2 * sigma2;
                epsilon = std::sqrt(eps[indexI] * eps[indexJ]); // TODO this can be cached
                d0 = x[indexI*3 + 0] - x[indexJ*3 + 0];
                d1 = x[indexI*3 + 1] - x[indexJ*3 + 1];
                d2 = x[indexI*3 + 2] - x[indexJ*3 + 2];
                dsqr = d0*d0 + d1*d1 + d2*d2;
                //check if is membrane -> need to skip attractive forces
                if (t[indexI] & 0x80000000 || t[indexJ] & 0x80000000) {
                    if (dsqr >= rt3_2 * sigma2) return;
                }

                l2NInvSquare = 1 / (dsqr);
                fac0 = 24 * epsilon * l2NInvSquare;
                l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;
                fac1_sum1 = sigma6 * l2NInvPow6;
                fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);

                auto access = _ff.access();
                access[indexI*3 + 0] -= fac0 * fac1 * d0;
                access[indexI*3 + 1] -= fac0 * fac1 * d1;
                access[indexI*3 + 2] -= fac0 * fac1 * d2;
                access[indexJ*3 + 0] += fac0 * fac1 * d0;
                access[indexJ*3 + 1] += fac0 * fac1 * d1;
                access[indexJ*3 + 2] += fac0 * fac1 * d2;
            });
            Kokkos::fence();
        });
    }
#endif
#endif
#endif
} // force