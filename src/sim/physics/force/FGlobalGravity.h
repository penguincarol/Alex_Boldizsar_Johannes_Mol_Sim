//
// Created by alex on 06.12.2022.
//

#pragma once

#include "ForceFunctorBase.h"

namespace sim::physics::force {
    /**
     * applies gGrav on all particles
     * */
    class FGlobalGravity : public ForceFunctorBase {
    private:
        double gGrav0, gGrav1, gGrav2;

    public:
        /**
         * the created instance will take ownership of ff and will delete it upon deconstruction.
         * */
        FGlobalGravity(double st,
                       double et,
                       double dt,
                       double eps,
                       double sig,
                       double gG0, double gG1, double gG2,
                       ParticleContainer &pc
        ) : ForceFunctorBase(st, et, dt, eps, sig, pc), gGrav0(gG0),gGrav1(gG1),gGrav2(gG2) {}

        void operator()() override;

        void setParticleContainer(ParticleContainer& pc) override;

        pair_fun_t& getForceFunction() override;

        fpair_fun_t getFastForceFunction() override;

        fpair_fun_alt_t getFastForceAltFunction() override;
    };

} // sim::physics::force
