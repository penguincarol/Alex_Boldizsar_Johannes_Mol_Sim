//
// Created by alex on 11.01.23.
//

#pragma once

#include "ForceFunctorBase.h"

namespace sim::physics::force {

    /**
     * Simulates the interaction of particles within a membrane (i.e. spring forces)
     * */
    class FMembrane : public ForceFunctorBase {
    public:
        /**
         * the created instance will take ownership of ff and will delete it upon deconstruction.
         * */
        FMembrane(double st,
                       double et,
                       double dt,
                       double eps,
                       double sig,
                       ParticleContainer &pc
        ) : ForceFunctorBase(st, et, dt, eps, sig, pc) {}

        void operator()() override;

        void setParticleContainer(ParticleContainer& pc) override;

        pair_fun_t& getForceFunction() override;

        fpair_fun_t getFastForceFunction() override;
    };

} // force
