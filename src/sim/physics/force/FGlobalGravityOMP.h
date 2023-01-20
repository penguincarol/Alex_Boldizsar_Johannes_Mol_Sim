//
// Created by alex on 11.01.23.
//

#pragma once

#include "ForceFunctorBase.h"

namespace sim::physics::force {

    /**
     * applies gGrav on all particles using OMP
     * */
    class FGlobalGravityOMP : public ForceFunctorBase {
    private:
        double gGrav0, gGrav1, gGrav2;

    public:
        /**
         * the created instance will take ownership of ff and will delete it upon deconstruction.
         * */
        FGlobalGravityOMP(double st,
                       double et,
                       double dt,
                       double eps,
                       double sig,
                       double gG0, double gG1, double gG2,
                       ParticleContainer &pc
        ) : ForceFunctorBase(st, et, dt, eps, sig, pc), gGrav0(gG0),gGrav1(gG1),gGrav2(gG2) {}

        void operator()() override;

        void setParticleContainer(ParticleContainer& pc) override;

    };

} // force

