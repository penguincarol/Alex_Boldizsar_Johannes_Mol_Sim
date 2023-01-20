//
// Created by alex on 11.01.23.
//

#pragma once

#include "ForceFunctorBase.h"
#include "FLennardJones.h"

namespace sim::physics::force {

    class FLennardJonesCellsOMP : public ForceFunctorBase {
    private:
        pair_fun_t pairFun;
        fpair_fun_t fpairFun;
        fpair_fun_alt_t fpairFunAlt;
        fpair_fun_ret_t fpairFunRet;
        FLennardJones forceDelegate;

        void setPairFun();

    public:
        fpair_fun_alt_t getFastForceAltFunction() override;

    public:
        FLennardJonesCellsOMP(double st,
                           double et,
                           double dt,
                           double eps,
                           double sig,
                           ParticleContainer &pc
        ) : ForceFunctorBase(st, et, dt, eps, sig, pc), forceDelegate(FLennardJones(st, et, dt, eps, sig, pc)) {
            setPairFun();
        }

        void operator()() override;

        void setParticleContainer(ParticleContainer& pc) override;

        pair_fun_t& getForceFunction() override;

        fpair_fun_t getFastForceFunction() override;

        fpair_fun_ret_t getFastForceRetFunction() override;
    };

} // force

