//
// Created by alex on 28.11.2022.
//

#pragma once

#include "BoundsFunctorBase.h"

namespace sim::physics::bounds {
    /**
     * Creates an outflowing bound on the side S
     * */
    template<sim::physics::bounds::side S>
    class BoundsOutflow : public BoundsFunctorBase<S> {
    public:
        ~BoundsOutflow() override = default;

        /**
         * @param st start time
         * @param et end time
         * @param dt delta time
         * @param eps epsilon
         * @param sig sigma
         * @param pc particle container
         * @param eOMP enable OMP flag
         * */
        BoundsOutflow(double st, double et, double dt, double eps, double sig, ParticleContainer &pc, bool eOMP)
                : BoundsFunctorBase<S>(st, et, dt, eps, sig, pc, eOMP) {}

        /**Removes particles upon crossing the boundary.*/
        void operator()() override {
            std::unordered_set<unsigned long> indices;
            this->particleContainer.template popExternalParticles<S>(indices);
            this->particleContainer.deactivateParticles(indices);
        };
    };
} // sim::physics::bounds

