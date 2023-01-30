//
// Created by alex on 11.01.23.
//

#include "FLennardJonesCellsOMP.h"
#include "defaults.h"

#include <iostream>

/**
* I am not writing yet another task initializer so this method just creates the 1d tasks using the 2d tasks
*/
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &orig)
{
    std::vector<T> ret;
    for(const auto &v: orig)
        ret.insert(ret.end(), v.begin(), v.end());
    return ret;
}

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

    #ifndef ONE_DIMENSIONAL_TASKS
    #ifndef THREE_DIMENSIONAL_TASKS
    #ifndef TASK_ORIENTED_2D
    void FLennardJonesCellsOMP::operator()() {
        //io::output::loggers::general->error("Default 2d thread oriented is used");
        particleContainer.runOnDataCell([&](std::vector<double> &force,
                                            std::vector<double> &oldForce,
                                            std::vector<double> &x,
                                            std::vector<double> &v,
                                            std::vector<double> &m,
                                            std::vector<int> &type,
                                            unsigned long count,
                                            ParticleContainer::VectorCoordWrapper &cells,
                                            std::vector<double> &eps,
                                            std::vector<double> &sig) {
            auto fpairFun = this->fpairFun;
#pragma omp parallel for default(none) shared(cells, x, eps, sig, m, type, force, fpairFun) //reduction(+:interactions)
            for (size_t cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
                auto &cell = cells[cellIndex];
                for (size_t i = 0; i < cell.size(); i++) {
                    for (size_t j = i + 1; j < cell.size(); j++) {
                        fpairFun(force, x, eps, sig, m, type, cell[i], cell[j]);
                    }
                }
            }
        });

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
            static const double rt3_2 = std::pow(2, 1.0 / 3.0);
            const std::vector<std::vector<std::pair<unsigned long, unsigned long>>> &alternativeTaskGroups = particleContainer.generateDistinctAlternativeCellNeighbours();
            double *_force = force.data();
            size_t size = force.size();
            double sigma, sigma2, sigma6, epsilon, dsqr, d0, d1, d2, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1;
            unsigned long indexI;
            unsigned long indexJ;
            unsigned long indexII;
            unsigned long indexJJ;
            unsigned long indexX;
            unsigned long indexY;
            unsigned long indexC0;
            unsigned long indexC1;

            size_t maxThreads = omp_get_max_threads();

#pragma omp parallel \
                default(none) \
                shared(size, x, t, cells, eps, sig, alternativeTaskGroups, fpairFun, force, m, maxThreads) \
                private(sigma, sigma2, sigma6, epsilon, dsqr, l2NInvSquare, d0, d1, d2, \
                        fac0, l2NInvPow6, fac1_sum1, fac1, indexI, indexJ, indexII, indexJJ, indexY, indexC0, indexC1) \
                firstprivate(indexX, rt3_2) \
                reduction(+:_force[:size])
            {
                //generate tasks: for all distinct cell neighbours
#pragma omp for simd
                for (indexX = 0; indexX < maxThreads; indexX++) {
                    for (indexY = 0; indexY < alternativeTaskGroups[indexX].size(); indexY++) {
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
    #endif
    #endif
    #endif
} // force