//
// Created by alex on 11.01.23.
//

#include "FLennardJonesCellsOMP.h"
#include "defaults.h"

#include <iostream>

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
    }

    fpair_fun_t FLennardJonesCellsOMP::getFastForceFunction() {
        return fpairFun;
    }

    void FLennardJonesCellsOMP::operator()() {
        particleContainer.runOnDataCell([&](std::vector<double> &force,
                                            std::vector<double> &oldForce,
                                            std::vector<double> &x,
                                            std::vector<double> &v,
                                            std::vector<double> &m,
                                            std::vector<int> &t,
                                            unsigned long count,
                                            ParticleContainer::VectorCoordWrapper &cells,
                                            std::vector<double> &eps,
                                            std::vector<double> &sig) {
            using cell_ptr = std::vector<unsigned long> *;
            std::vector<std::vector<cell_ptr>> tasks;
            auto fpairFun = this->fpairFun;
            cell_ptr c_ptr;
            static const double rt3_2 = std::pow(2, 1 / 3);
            const std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>>& taskGroups = particleContainer.generateDistinctCellNeighbours();
            const std::vector<std::vector<std::pair<unsigned long, unsigned long>>>& alternativeTaskGroups = particleContainer.generateDistinctAlternativeCellNeighbours();
            double *f = force.data();
            size_t size = force.size();
            double sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1;
            unsigned long indexI;
            unsigned long indexJ;
            unsigned long indexII;
            unsigned long indexJJ;
            unsigned long indexX;
            unsigned long indexY;
            unsigned long indexZ;
            unsigned long indexC0;
            unsigned long indexC1;

            size_t maxThreads = omp_get_max_threads();
            std::vector<size_t> interactions;
            interactions.clear();
            tasks.clear();
            interactions.resize(maxThreads,0);
            tasks.resize(maxThreads);

            for (auto &cell: cells) {
                size_t indexMin = 0;
                size_t last_count = interactions[indexMin];
                for(size_t i = 1; i < maxThreads; i++) {
                    if(interactions[i] <= last_count) {
                        last_count = interactions[i];
                        indexMin = i;
                    }
                }

                tasks[indexMin].emplace_back(&cell);
                interactions[indexMin] += cell.size()*(cell.size()-1)/2;
            }

//            for(size_t i= 0; i < maxThreads; i++) {
//                std::cout << "Task " << i << ": " << interactions[i] << std::endl;
//            }
            #pragma omp parallel default(none) shared(force, x, m, t, eps, sig, tasks, fpairFun) private(c_ptr, indexI, indexJ,indexII,indexJJ,indexX,indexY)
            //#pragma omp parallel default(none) shared(size, force, x, v, m, t, count, cells, eps, sig, tasks, buffer_size, fpairFun, rt3_2, taskGroups) private(c_ptr, sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1, indexI, indexJ,indexII,indexJJ,indexX,indexY,indexC0,indexC1,indexZ)
            {
                #pragma omp for
                for (indexII = 0; indexII < omp_get_max_threads(); indexII++) {
                    for (indexJJ = 0; indexJJ < tasks[indexII].size(); indexJJ++) {
                        c_ptr = tasks[indexII][indexJJ];
                        for (indexX = 0; indexX < c_ptr->size(); indexX++) {
                            for (indexY = indexX + 1; indexY < c_ptr->size(); indexY++) {
                                indexI = (*c_ptr)[indexX];
                                indexJ = (*c_ptr)[indexY];
                                fpairFun(force, x, eps, sig, m, t, indexI, indexJ);
                            }
                        }
                    }
                }

                #pragma omp barrier
            }

//            for(size_t i= 0; i < maxThreads; i++) {
//                std::cout << "Task " << i << ": " << alternativeTaskGroups[i].size() << std::endl;
//            }

            #pragma omp parallel \
                default(none) \
                shared(size, x, t, cells, eps, sig,alternativeTaskGroups,fpairFun,force,m) \
                private(sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, \
                        fac0, l2NInvPow6, fac1_sum1, fac1, indexI, indexJ,indexII,indexJJ,indexY,indexC0,indexC1) \
                firstprivate(indexX, rt3_2) \
                reduction(+:f[:size])
            {
                //generate tasks: for all distinct cell neighbours
                #pragma omp for
                for(indexX = 0; indexX < omp_get_max_threads(); indexX++) {
                    for(indexY = 0; indexY < alternativeTaskGroups[indexX].size(); indexY++){
                        indexC0 = alternativeTaskGroups[indexX][indexY].first;
                        indexC1 = alternativeTaskGroups[indexX][indexY].second;
                        for (indexII = 0; indexII < cells[indexC0].size(); indexII++) {
                            for (indexJJ = 0; indexJJ < cells[indexC1].size(); indexJJ++) {
                                indexI = cells[indexC0][indexII];
                                indexJ = cells[indexC1][indexJJ];
                                //fpairFun(force, x, eps, sig, m, t, indexI, indexJ);
                                sigma = (sig[indexI] + sig[indexJ]) / 2;
                                sigma2 = sigma * sigma;
                                sigma6 = sigma2 * sigma2 * sigma2;
                                epsilon = std::sqrt(eps[indexI] * eps[indexJ]); // TODO this can be cached
                                d0 = x[indexI * 3 + 0] - x[indexJ * 3 + 0];
                                d1 = x[indexI * 3 + 1] - x[indexJ * 3 + 1];
                                d2 = x[indexI * 3 + 2] - x[indexJ * 3 + 2];
                                dsqr = d0 * d0 + d1 * d1 + d2 * d2;
                                //check if is membrane -> need to skip attractive forces
                                if (t[indexI] & 0x80000000 || t[indexJ] & 0x80000000) {
                                    if (dsqr >= rt3_2 * sigma2) continue;
                                }

                                l2NInvSquare = 1 / (dsqr);
                                fac0 = 24 * epsilon * l2NInvSquare;
                                l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;
                                fac1_sum1 = sigma6 * l2NInvPow6;
                                fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);

                                f[indexI * 3 + 0] -= fac0 * fac1 * d0;
                                f[indexI * 3 + 1] -= fac0 * fac1 * d1;
                                f[indexI * 3 + 2] -= fac0 * fac1 * d2;
                                f[indexJ * 3 + 0] += fac0 * fac1 * d0;
                                f[indexJ * 3 + 1] += fac0 * fac1 * d1;
                                f[indexJ * 3 + 2] += fac0 * fac1 * d2;
                            }
                        }
                    }

                }
            }
        });
    }
} // force