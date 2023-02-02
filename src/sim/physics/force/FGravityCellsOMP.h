//
// Created by alex on 11.01.23.
//

#pragma once

#include "ForceFunctorBase.h"
#include "FGravity.h"

namespace sim::physics::force {
    /**
     * Performs gravity interaction between two particles while using OMP and the linked cell approach.
     * WARNING: this is not supported any more
     * */
    class FGravityCellsOMP : public ForceFunctorBase{
    private:
        pair_fun_t pairFun;
        fpair_fun_t fpairFun;
        fpair_fun_alt_t fpairFunAlt;
        fpair_fun_ret_t fpairFunRet;
        FGravity forceDelegate;

        void setPairFun();

    public:
        fpair_fun_ret_t getFastForceRetFunction() override;

    public:
        /**
         * @param st start time
         * @param et end time
         * @param dt delta time
         * @param eps epsilon
         * @param sig sigma
         * @param pc particle container
         * */
        FGravityCellsOMP(double st,
                      double et,
                      double dt,
                      double eps,
                      double sig,
                      ParticleContainer &pc
        ) : ForceFunctorBase(st, et, dt, eps, sig, pc), forceDelegate(st, et, dt, eps, sig, pc) {
            std::cerr << "Sorry not supported, please do not use Gravity with Cells and OMP" << std::endl;
            exit(-1);
            //setPairFun();
        }

        void operator()() override;

        void setParticleContainer(ParticleContainer& pc) override;

        pair_fun_t& getForceFunction() override;

        fpair_fun_t getFastForceFunction() override;

        fpair_fun_alt_t getFastForceAltFunction() override;
    };

} // force

