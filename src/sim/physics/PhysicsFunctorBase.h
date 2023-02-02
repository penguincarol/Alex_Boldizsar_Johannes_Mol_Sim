//
// Created by alex on 26.11.2022.
//

#pragma once

#include "data/ParticleContainer.h"

namespace sim::physics {
    /**
     * Base class for all functors.
     * Provides storage for basic information about the current simulation.
     * */
    class PhysicsFunctorBase {
    protected:
        const double start_time;
        const double end_time;
        const double delta_t;
        const double epsilon;
        const double sigma;
        ParticleContainer &particleContainer;

    public:
        PhysicsFunctorBase() = delete;

        /**
         * @param st start time
         * @param et end time
         * @param dt delta time
         * @param eps epsilon
         * @param sig sigma
         * @param pc particle container
         * */
        PhysicsFunctorBase(double st, double et, double dt, double eps, double sig, ParticleContainer& pc) :
                start_time(st), end_time(et), delta_t(dt), epsilon(eps), sigma(sig), particleContainer(pc) {}

        virtual ~PhysicsFunctorBase() = default;

        virtual void operator()() = 0;

        virtual void setParticleContainer(ParticleContainer& pc) { particleContainer = pc; }
    };

} // sim::physics