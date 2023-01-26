//
// Created by alex on 11.01.23.
//

#include "FGravityCellsOMP.h"
#include "defaults.h"

namespace sim::physics::force {
    void FGravityCellsOMP::setPairFun() {
        pairFun = forceDelegate.getForceFunction();
        fpairFun = forceDelegate.getFastForceFunction();
        fpairFunAlt = forceDelegate.getFastForceAltFunction();
        fpairFunRet = forceDelegate.getFastForceRetFunction();
    }

    void FGravityCellsOMP::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
        forceDelegate.setParticleContainer(pc);
        setPairFun();
    }

    pair_fun_t &FGravityCellsOMP::getForceFunction() {
        return pairFun;
    }

    fpair_fun_t FGravityCellsOMP::getFastForceFunction() {
        return fpairFun;
    }

    void FGravityCellsOMP::operator()() {
        particleContainer.runOnDataCell([&](Kokkos::View<double*> &force,
                                            Kokkos::View<double*> &oldForce,
                                            Kokkos::View<double*> &x,
                                            Kokkos::View<double*> &v,
                                            Kokkos::View<double*> &m,
                                            Kokkos::View<int*> &type,
                                            unsigned long count,
                                            ParticleContainer::VectorCoordWrapper& cells,
                                            Kokkos::View<double*> &eps,
                                            Kokkos::View<double*> &sig){
            // TODO impl
        });
    }

    fpair_fun_alt_t FGravityCellsOMP::getFastForceAltFunction() {
        return fpairFunAlt;
    }

    fpair_fun_ret_t FGravityCellsOMP::getFastForceRetFunction() {
        return fpairFunRet;
    }
} // force