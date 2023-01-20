//
// Created by alex on 29.11.2022.
//

#pragma once

#include <functional>

#include "sim/physics/PhysicsFunctorBase.h"

namespace sim::physics::force {
    using pair_fun_t = std::function<void(Particle &, Particle &)>;
    using fpair_fun_t = void(*)(std::vector<double> &force,
                                   std::vector<double> &x,
                                   std::vector<double> &eps,
                                   std::vector<double> &sig,
                                   std::vector<double> &m,
                                   std::vector<int> &t,
                                   unsigned long indexI, unsigned long indexJ);
    using fpair_fun_alt_t = void(*)(std::vector<double> &force,
                                    std::vector<double> &x,
                                    std::vector<double> &eps,
                                    std::vector<double> &sig,
                                    std::vector<double> &m,
                                    std::vector<int> &t,
                                    unsigned long indexI,
                                    double xJ0, double xJ1, double xJ2,
                                    double epsJ, double sigJ, double mJ, int tJ);
    using fpair_fun_ret_t = std::array<double,3>(*)(std::vector<double> &force,
                                                    std::vector<double> &x,
                                                    std::vector<double> &eps,
                                                    std::vector<double> &sig,
                                                    std::vector<double> &m,
                                                    std::vector<int> &t,
                                                    unsigned long indexI,
                                                    double xJ0, double xJ1, double xJ2,
                                                    double epsJ, double sigJ, double mJ, int tJ);

    class ForceFunctorBase : public PhysicsFunctorBase {
    public:
        ForceFunctorBase(double st, double et, double dt, double eps, double sig, ParticleContainer &pc)
                : PhysicsFunctorBase(st, et, dt, eps, sig, pc) {};

        /**
         * Returns the function that calculates the force between two particles
         * */
        virtual pair_fun_t& getForceFunction() {
            throw std::runtime_error{"not implemented"};
        }

        /**
         * Returns the fast version of the pairwise force function
         * */
        virtual fpair_fun_t getFastForceFunction() {
            throw std::runtime_error{"not implemented"};
        }

        /**
         * Returns the fast version of the pairwise force function, in which the second particle is a temporary particle
         * */
         virtual fpair_fun_alt_t getFastForceAltFunction() {
            throw std::runtime_error{"not implemented"};
        }

         /**
          * Same as FFAlt, but does not write force back, instead returns the value
          * */
         virtual fpair_fun_ret_t getFastForceRetFunction() {
             throw std::runtime_error{"not implemented"};
         }
    };
} // sim::physics::force


