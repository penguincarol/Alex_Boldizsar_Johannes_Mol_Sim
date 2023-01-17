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
            size_t buffer_size = 0;
            auto fpairFun = this->fpairFun;
            cell_ptr c_ptr;
            static const double rt3_2 = std::pow(2, 1 / 3);
            const std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>>& taskGroups = particleContainer.generateDistinctCellNeighbours();
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

            tasks.emplace_back();
            for (auto &cell: cells) {
                if (buffer_size > max_thread_tasks) {
                    tasks.emplace_back();
                    buffer_size = 0;
                }
                tasks.back().emplace_back(&cell);
                buffer_size += cell.size();
            }

            #pragma omp parallel default(none) shared(force, x, m, t, eps, sig, tasks, fpairFun) private(c_ptr, indexI, indexJ,indexII,indexJJ,indexX,indexY)
            //#pragma omp parallel default(none) shared(size, force, x, v, m, t, count, cells, eps, sig, tasks, buffer_size, fpairFun, rt3_2, taskGroups) private(c_ptr, sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1, indexI, indexJ,indexII,indexJJ,indexX,indexY,indexC0,indexC1,indexZ)
            {
                #pragma omp for
                for (indexII = 0; indexII < tasks.size(); indexII++) {
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

            for(indexX = 0; indexX < 26; indexX++) { //task group index: indexX
                #pragma omp parallel default(none) shared(force, x, m, t, cells, eps, sig, fpairFun, taskGroups) private(indexI, indexJ,indexII,indexJJ,indexX,indexY,indexC0,indexC1,indexZ)
                {
                    #pragma omp for schedule(static,1)
                    for (indexY = 0; indexY < taskGroups[indexX].size(); indexY++) { //task index: indexY - each thread gets one task
                        for (indexZ = 0; indexZ < taskGroups[indexX][indexY].size(); indexZ++) { //pair index: indexII
                            indexC0 = taskGroups[indexX][indexY][indexZ].first;
                            indexC1 = taskGroups[indexX][indexY][indexZ].second;
                            for (indexII = 0; indexII < cells[indexC0].size(); indexII++) {
                                for (indexJJ = 0; indexJJ < cells[indexC1].size(); indexJJ++) {
                                    indexI = cells[indexC0][indexII];
                                    indexJ = cells[indexC1][indexJJ];
                                    fpairFun(force, x, eps, sig, m, t, indexI, indexJ);
                                }
                            }
                        }
                    }
                    #pragma omp barrier
                }
            }

//            #pragma omp parallel default(none) shared(size, force, x, v, m, t, count, cells, eps, sig, tasks, buffer_size, fpairFun, rt3_2, taskGroups) private(c_ptr, sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1, indexI, indexJ,indexII,indexJJ,indexX,indexY,indexC0,indexC1,indexZ) reduction(+:f[:size])
//            {
//                //generate tasks: for all distinct cell neighbours
//                #pragma omp for
//                for(indexX = 0; indexX < alternativeTaskGroups.size(); indexX++) {
//                    indexC0 = alternativeTaskGroups[indexX].first;
//                    indexC1 = alternativeTaskGroups[indexX].second;
//                    for (indexII = 0; indexII < cells[indexC0].size(); indexII++) {
//                        for (indexJJ = 0; indexJJ < cells[indexC1].size(); indexJJ++) {
//                            indexI = cells[indexC0][indexII];
//                            indexJ = cells[indexC1][indexJJ];
//                            sigma = (sig[indexI] + sig[indexJ]) / 2;
//                            sigma2 = sigma * sigma;
//                            sigma6 = sigma2 * sigma2 * sigma2;
//                            epsilon = std::sqrt(eps[indexI] * eps[indexJ]); // TODO this can be cached
//                            d0 = x[indexI * 3 + 0] - x[indexJ * 3 + 0];
//                            d1 = x[indexI * 3 + 1] - x[indexJ * 3 + 1];
//                            d2 = x[indexI * 3 + 2] - x[indexJ * 3 + 2];
//                            dsqr = d0 * d0 + d1 * d1 + d2 * d2;
//                            //check if is membrane -> need to skip attractive forces
//                            if (t[indexI] & 0x80000000 || t[indexJ] & 0x80000000) {
//                                if (dsqr >= rt3_2 * sigma2) continue;
//                            }
//
//                            l2NInvSquare = 1 / (dsqr);
//                            fac0 = 24 * epsilon * l2NInvSquare;
//                            l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;
//                            fac1_sum1 = sigma6 * l2NInvPow6;
//                            fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);
//
//                            f[indexI * 3 + 0] -= fac0 * fac1 * d0;
//                            f[indexI * 3 + 1] -= fac0 * fac1 * d1;
//                            f[indexI * 3 + 2] -= fac0 * fac1 * d2;
//                            f[indexJ * 3 + 0] += fac0 * fac1 * d0;
//                            f[indexJ * 3 + 1] += fac0 * fac1 * d1;
//                            f[indexJ * 3 + 2] += fac0 * fac1 * d2;
//                        }
//                    }
//
//                }
//            }
        });
    }
} // force