//
// Created by alex on 26.11.2022.
//

#pragma once

#include "ForceFunctorBase.h"

namespace sim::physics::force {
    /**
     * calculate the force for all particles using the Lennard-Jones potential
     * */
    class FLennardJones : public ForceFunctorBase {
    private:
        pair_fun_t pairFun;
        fpair_fun_t fpairFun;

        void setPairFun();

    public:
        /**
         * @param st start time
         * @param et end time
         * @param dt delta time
         * @param eps epsilon
         * @param sig sigma
         * @param pc particle container
         * */
        FLennardJones(double st,
                      double et,
                      double dt,
                      double eps,
                      double sig,
                      ParticleContainer &pc
        ) : ForceFunctorBase(st, et, dt, eps, sig, pc) {
            setPairFun();
        }

        void operator()() override;

        void setParticleContainer(ParticleContainer& pc) override;

        pair_fun_t& getForceFunction() override;

        fpair_fun_ret_t getFastForceRetFunction() override;

        fpair_fun_alt_t getFastForceAltFunction() override;

        fpair_fun_t getFastForceFunction() override;
    };

} // sim::physics::force