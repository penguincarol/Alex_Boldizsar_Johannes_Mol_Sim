//
// Created by alex on 11.01.23.
//

#pragma once

#include "ForceFunctorBase.h"

namespace sim::physics::force {
    /**
     * Simulates a pulling force on the specified coordinates within a membrane
     * */
    class FMembranePull : public ForceFunctorBase{
    private:
        double current_time;
    public:
        /**
         * @param st start time
         * @param et end time
         * @param dt delta time
         * @param eps epsilon
         * @param sig sigma
         * @param pc particle container
         * */
        FMembranePull(double st,
                     double et,
                     double dt,
                     double eps,
                     double sig,
                     ParticleContainer &pc
        ) : ForceFunctorBase(st, et, dt, eps, sig, pc), current_time(st) {}

        void operator()() override;

        void setParticleContainer(ParticleContainer& pc) override;

    };

} // force
