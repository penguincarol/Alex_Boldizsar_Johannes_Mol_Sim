//
// Created by alex on 11.01.23.
//

#pragma once

#include "ForceFunctorBase.h"
#include "FGravity.h"

namespace sim::physics::force {

    /**
     * Performs gravity interaction between two particles while using OMP.
     * */
    class FGravityOMP : public ForceFunctorBase {
    private:
        pair_fun_t pairFun;
        fpair_fun_t fpairFun;
        FGravity forceDelegate;

        void setPairFun();

    public:
        FGravityOMP(double st,
                      double et,
                      double dt,
                      double eps,
                      double sig,
                      ParticleContainer &pc
        ) : ForceFunctorBase(st, et, dt, eps, sig, pc), forceDelegate(st, et, dt, eps, sig, pc) {
            setPairFun();
        }

        void operator()() override;

        void setParticleContainer(ParticleContainer& pc) override;

        pair_fun_t& getForceFunction() override;

        fpair_fun_t getFastForceFunction() override;
    };

} // force

