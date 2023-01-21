//
// Created by alex on 11.01.23.
//

#include "FLennardJonesCellsOMP.h"
#include "defaults.h"

#include <iostream>
#include <immintrin.h>

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
            //const std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>>& taskGroups = particleContainer.generateDistinctCellNeighbours();
            const std::vector<std::vector<std::pair<unsigned long, unsigned long>>>& alternativeTaskGroups = particleContainer.generateDistinctAlternativeCellNeighbours();
            double *_force = force.data();
            double *_x = x.data();
            double *_sig = sig.data();
            double *_eps = eps.data();
            size_t size = force.size();
            double sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1;
            unsigned long indexI;
            unsigned long indexJ;
            unsigned long indexII;
            unsigned long indexJJ;
            unsigned long indexX;
            unsigned long indexY;
            //unsigned long indexZ;
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
            #pragma omp parallel default(none) shared(force, x, m, t, eps, sig, tasks, fpairFun, maxThreads) private(c_ptr, indexI, indexJ,indexII,indexJJ,indexX,indexY)
            //#pragma omp parallel default(none) shared(size, force, x, v, m, t, count, cells, eps, sig, tasks, buffer_size, fpairFun, rt3_2, taskGroups) private(c_ptr, sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1, indexI, indexJ,indexII,indexJJ,indexX,indexY,indexC0,indexC1,indexZ)
            {
                #pragma omp for
                for (indexII = 0; indexII < maxThreads; indexII++) {
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
                shared(size, _x, x, t, cells, _eps, eps, _sig, sig,alternativeTaskGroups,fpairFun,force,m, maxThreads) \
                private(sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, \
                        fac0, l2NInvPow6, fac1_sum1, fac1, indexI, indexJ,indexII,indexJJ,indexY,indexC0,indexC1) \
                firstprivate(indexX, rt3_2) \
                reduction(+:_force[:size])
            {
                //generate tasks: for all distinct cell neighbours
                #pragma omp for
                for(indexX = 0; indexX < maxThreads; indexX++) {
                    for(indexY = 0; indexY < alternativeTaskGroups[indexX].size(); indexY++){
                        indexC0 = alternativeTaskGroups[indexX][indexY].first;
                        indexC1 = alternativeTaskGroups[indexX][indexY].second;
                        //TASK:
                        //first do 0-3 iterations for cellI s.t. its size%4==0
                        //handle cellJ: do make size%4==0, then vectorize
                        //then handle rest of cellI vectorized
                        //handle cellJ in the same way
                        /*
                         * static const __m256i xMask = _mm256_set_epi64x(0, -1, -1, -1);
                                __m256d xI = _mm256_maskload_pd(_x + 2*indexI + indexI, xMask);
                                __m256d xJ = _mm256_maskload_pd(_x + 2*indexJ + indexJ, xMask);
                                __m256d d = _mm256_sub_pd(xI, xJ);
                                //for sigma/epsilon to permute: _mm256_permute_pd -> for 4 particles once load all values
                         * */

                        //init
                        indexII = 0;
                        indexJJ = 0;
                        static const __m256i xMask = _mm256_set_epi64x(0, -1, -1, -1);
                        static const __m256d half = _mm256_set1_pd(0.5);

                        // make cellI size divisible by 4
                        for(; indexII < cells[indexC0].size() % 4; indexII++) {
                            indexI = cells[indexC0][indexII];
                            double sigIs = sig[indexI];
                            double epsIs = eps[indexI];
                            __m256d xI = _mm256_maskload_pd(_x + 2 * indexI + indexI, xMask);
                            __m256d sigI = _mm256_set1_pd(sigIs); //sigma of I in all positions
                            __m256d epsI = _mm256_set1_pd(epsIs); //epsilon of I in all positions

                            //I is not vector, J is not vector
                            for(; indexJJ < cells[indexC1].size() % 4; indexJJ++) {
                                indexJ = cells[indexC1][indexJJ];
                                sigma = (sigIs + sig[indexJ]) / 2;
                                epsilon = std::sqrt(epsIs * eps[indexJ]);
                                sigma6 = std::pow(sigma,6);

                                __m256d xJ = _mm256_maskload_pd(_x + 2 * indexJ + indexJ, xMask);
                                __m256d d = _mm256_sub_pd(xI, xJ);
                                __m256d d_sqr = _mm256_mul_pd(d,d);
                                __m256d hadd = _mm256_hadd_pd(d_sqr, d_sqr); //dq0+dq1 -> [63:0], dq2+0 -> [191:128]
                                __m128d upper = _mm256_extractf128_pd(hadd, 1);
                                __m128d d_sum = _mm_add_pd(upper, _mm256_castpd256_pd128(hadd));
                                dsqr = _mm_cvtsd_f64(d_sum);
                                if (t[indexI] & 0x80000000 || t[indexJ] & 0x80000000) {
                                    if (dsqr >= rt3_2 * sigma2) continue;
                                }

                                l2NInvSquare = 1 / (dsqr);
                                fac0 = 24 * epsilon * l2NInvSquare;
                                l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;
                                fac1_sum1 = sigma6 * l2NInvPow6;
                                fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);
                                _force[indexI * 3 + 0] -= fac0 * fac1 * d0;
                                _force[indexI * 3 + 1] -= fac0 * fac1 * d1;
                                _force[indexI * 3 + 2] -= fac0 * fac1 * d2;
                                _force[indexJ * 3 + 0] += fac0 * fac1 * d0;
                                _force[indexJ * 3 + 1] += fac0 * fac1 * d1;
                                _force[indexJ * 3 + 2] += fac0 * fac1 * d2;

                            }
                            //I is not vector, J is vector now
                            for(; indexJJ < cells[indexC1].size(); indexJJ += 4) {
                                indexJ = cells[indexC1][indexJJ];
                                __m256d sigJ = _mm256_loadu_pd(_sig + indexJ); // sigma of all 4 particles
                                __m256d epsJ = _mm256_loadu_pd(_eps + indexJ); // epsilon of all 4 particles
                                __m256d xJ0 = _mm256_maskload_pd(_x + 2 * indexJ + indexJ, xMask);
                                __m256d xJ1 = _mm256_maskload_pd(_x + 2 * indexJ + indexJ + 2 + 1, xMask);
                                __m256d xJ2 = _mm256_maskload_pd(_x + 2 * indexJ + indexJ + 4 + 2, xMask);
                                __m256d xJ3 = _mm256_maskload_pd(_x + 2 * indexJ + indexJ + 8 + 1, xMask);

                                __m256d tmpSig = _mm256_add_pd(sigI, sigJ);
                                __m256d sig = _mm256_mul_pd(tmpSig, half);
                                __m256d tmpEps = _mm256_mul_pd(epsI, epsJ);
                                __m256d eps = _mm256_sqrt_pd(tmpEps);

                                __m256d d0 = _mm256_sub_pd(xI, xJ0);
                                __m256d d1 = _mm256_sub_pd(xI, xJ1);
                                __m256d d2 = _mm256_sub_pd(xI, xJ2);
                                __m256d d3 = _mm256_sub_pd(xI, xJ3);

                                __m256d d = _mm256_sub_pd(xI, xJ);
                                __m256d d_sqr = _mm256_mul_pd(d,d);
                                __m256d hadd = _mm256_hadd_pd(d_sqr, d_sqr); //dq0+dq1 -> [63:0], dq2+0 -> [191:128]
                                __m128d upper = _mm256_extractf128_pd(hadd, 1);
                                __m128d d_sum = _mm_add_pd(upper, _mm256_castpd256_pd128(hadd));
                            }
                        }
                        // cellI is divisible by 4 now
                        for(; indexII < cells[indexC0].size(); indexII += 4){
                            //I is vector, J is not necessarily vector
                            for(; indexJJ < cells[indexC1].size() % 4; indexJJ++) {
                                //do calc
                            }
                            //I not vector, J is vector
                            for(; indexJJ < cells[indexC1].size(); indexJJ += 4) {
                                //do calc
                            }
                        }

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

                                _force[indexI * 3 + 0] -= fac0 * fac1 * d0;
                                _force[indexI * 3 + 1] -= fac0 * fac1 * d1;
                                _force[indexI * 3 + 2] -= fac0 * fac1 * d2;
                                _force[indexJ * 3 + 0] += fac0 * fac1 * d0;
                                _force[indexJ * 3 + 1] += fac0 * fac1 * d1;
                                _force[indexJ * 3 + 2] += fac0 * fac1 * d2;
                            }
                        }
                    }

                }
            }
        });
    }
} // force