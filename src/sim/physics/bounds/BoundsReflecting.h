//
// Created by alex on 28.11.2022.
//

#pragma once

#include "BoundsFunctorBase.h"
#include "sim/physics/force/ForceHandler.h"

namespace sim::physics::bounds {
    /**
     * Creates a reflecting bound on the side S
     * */
    template<sim::physics::bounds::side S>
    class BoundsReflecting : public BoundsFunctorBase<S> {
    private:
        sim::physics::force::ForceHandler& forceHandler;

    public:
        ~BoundsReflecting() override = default;

        /**
         * @brief Creates a reflecting bound with the specified parameters
         * @param st start time
         * @param et end time
         * @param dt delta time
         * @param eps epsilon
         * @param sig sigma
         * @param pc particle container
         * @param fh force handler
         * @param eOMP enable OMP flag
         */
        BoundsReflecting(double st, double et, double dt, double eps, double sig, ParticleContainer &pc, sim::physics::force::ForceHandler& fh, bool eOMP)
                : BoundsFunctorBase<S>(st, et, dt, eps, sig, pc, eOMP), forceHandler(fh) {}

        /**Reflects particle upon nearing the border.*/
        void operator()() override {
            this->particleContainer.template forEachParticleHaloPairInSide<S>(forceHandler.getFastForceAltFunction());
        }
    };
} // sim::physics::bounds
