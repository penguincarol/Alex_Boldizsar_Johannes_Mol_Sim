//
// Created by alex on 11.01.23.
//

#pragma once

#include "ForceFunctorBase.h"
#include "FLennardJones.h"

namespace sim::physics::force {

    class FLennardJonesCellsOMP : public ForceFunctorBase {
    private:
        pair_fun_t pairFun;
        fpair_fun_t fpairFun;
        fpair_fun_alt_t fpairFunAlt;
        fpair_fun_ret_t fpairFunRet;
        FLennardJones forceDelegate;

        void setPairFun();

    public:
        fpair_fun_alt_t getFastForceAltFunction() override;

    public:
        FLennardJonesCellsOMP(double st,
                           double et,
                           double dt,
                           double eps,
                           double sig,
                           ParticleContainer &pc
        ) : ForceFunctorBase(st, et, dt, eps, sig, pc), forceDelegate(FLennardJones(st, et, dt, eps, sig, pc)) {
            setPairFun();
        }

        /**
         * Calculates the Lennard Jones Forces using the Linked Cell algorithm.
         *
         * There are 4 implementations of this method using the following task models: 1D tasks, 2D color oriented, 2D thread oriented and 3D oriented.
         * For more information look into the presentation in the branch Jo/presentation
         *
         * The 2D thread oriented approach and the 3D approach can both utilize the round robin distribution strategy or the greedy strategy.
         *
         * By default the 2D thread oriented approach using the greedy distribution strategy gets chosen.
         *
         * The implementation used can be chosen by adding the flags '''-Dthree_dim_tasks=1''', '''-Done_dim_tasks=1''', '''-Dtask_oriented_2d=1''' and
         * '''-Dround_robin_distr=1''' when running CMake (e.g. cmake -Dthree_dim_tasks=1 -Dround_robin_distr=1 ..)
         *
         * \image html speedup_2D_st_comp.png width=1px
         * \image html speedup_1D_st_comp.png width=1px
         * \image html speedup_3D_st_comp.png width=1px
         *
         * \htmlonly
         * <div style="display: flex;">
         *     <div style="flex: 1;">
         *     <img src="speedup_1D_st_comp.png"/>
         *     </div>
         *     <div style="flex: 1;">
         *         <img src="speedup_2D_st_comp.png"/>
         *     </div>
         *      <div style="flex: 1;">
         *         <img src="speedup_1D_st_comp.png"/>
         *     </div>
         * </div>
         * \endhtmlonly
         *
         */
        void operator()() override;

        void setParticleContainer(ParticleContainer& pc) override;

        pair_fun_t& getForceFunction() override;

        fpair_fun_t getFastForceFunction() override;

        fpair_fun_ret_t getFastForceRetFunction() override;
    };

} // force

