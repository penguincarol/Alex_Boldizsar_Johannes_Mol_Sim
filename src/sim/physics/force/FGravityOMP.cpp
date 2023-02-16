//
// Created by alex on 11.01.23.
//

#include "FGravityOMP.h"

namespace sim::physics::force {
    void FGravityOMP::setPairFun() {
        pairFun = forceDelegate.getForceFunction();
        fpairFun = forceDelegate.getFastForceFunction();
        fpairFunAlt = forceDelegate.getFastForceAltFunction();
        fpairFunRet = forceDelegate.getFastForceRetFunction();
    }

    void FGravityOMP::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
        forceDelegate.setParticleContainer(pc);
        setPairFun();
    }

    pair_fun_t &FGravityOMP::getForceFunction() {
        return pairFun;
    }

    fpair_fun_t FGravityOMP::getFastForceFunction() {
        return fpairFun;
    }

    void FGravityOMP::operator()() {
        particleContainer.runOnActiveData([](Kokkos::View<double*> &force,
                                       Kokkos::View<double*> &oldForce,
                                       Kokkos::View<double*> &x,
                                       Kokkos::View<double*> &v,
                                       Kokkos::View<double*> &m,
                                       Kokkos::View<int*> &type,
                                       unsigned long count,
                                       Kokkos::View<double*> &eps,
                                       Kokkos::View<double*> &sig,
                                       std::vector<unsigned long> &activeParticles) {

            double d0, d1, d2, scalar;
            unsigned long indexI;
            unsigned long indexJ;
            unsigned long size = activeParticles.size();
            unsigned long endIndex = size * (size + 1) / 2;
            double* f = force.data();

#pragma omp parallel default(none) shared(force, oldForce, x, m, count, endIndex, f) private(d0, d1, d2, indexI, indexJ, scalar)
            {
#pragma omp for reduction(+:f[:count*3])
                for(unsigned long globalIndex = 0; globalIndex < endIndex; globalIndex++){
                    indexI = globalIndex / count;
                    indexJ = globalIndex % count;
                    if (indexI > indexJ) {
                        indexI = count - indexI;
                        indexJ = count - indexJ - 1;
                    }
                    if(indexI == indexJ) continue;
                    d0 = x[indexI*3 + 0] - x[indexJ*3 + 0];
                    d1 = x[indexI*3 + 1] - x[indexJ*3 + 1];
                    d2 = x[indexI*3 + 2] - x[indexJ*3 + 2];
                    scalar = m[indexI] * m[indexJ] * std::pow(1.0/std::sqrt(d0*d0+d1*d1+d2*d2),3);


                    f[indexI*3 + 0] -= d0 * scalar;
                    f[indexI*3 + 1] -= d1 * scalar;
                    f[indexI*3 + 2] -= d2 * scalar;
                    f[indexJ*3 + 0] += d0 * scalar;;
                    f[indexJ*3 + 1] += d1 * scalar;;
                    f[indexJ*3 + 2] += d2 * scalar;;
                }
            }
        });
    }

    fpair_fun_alt_t FGravityOMP::getFastForceAltFunction() {
        return fpairFunAlt;
    }

    fpair_fun_ret_t FGravityOMP::getFastForceRetFunction() {
        return fpairFunRet;
    }

} // force