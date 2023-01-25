//
// Created by alex on 26.11.2022.
//

#include "FGravity.h"

namespace sim::physics::force {
    void FGravity::operator()() {
        particleContainer.forAllPairs(pairFun);
    }

    pair_fun_t &FGravity::getForceFunction() {
        return pairFun;
    }

    static void fastPairFunction(std::vector<double> &force,
                                 std::vector<double> &x,
                                 std::vector<double> &eps,
                                 std::vector<double> &sig,
                                 std::vector<double> &m,
                                 std::vector<int> &t,
                                 unsigned long indexI, unsigned long indexJ) {
        double d0, d1, d2, s, df0, df1, df2;
        d0 = x[3*indexI + 0] - x[3*indexJ + 0];
        d1 = x[3*indexI + 1] - x[3*indexJ + 1];
        d2 = x[3*indexI + 2] - x[3*indexJ + 2];
        s = m[indexI] * m[indexJ] * std::pow(1.0 / std::sqrt(d0 * d0 + d1 * d1 + d2 * d2), 3);
        df0 = d0 * s;
        df1 = d1 * s;
        df2 = d2 * s;
        force[3*indexI + 0] -= df0;
        force[3*indexI + 1] -= df1;
        force[3*indexI + 2] -= df2;
        force[3*indexJ + 0] += df0;
        force[3*indexJ + 1] += df1;
        force[3*indexJ + 2] += df2;
    }

    void FGravity::setPairFun() {
        pairFun = [&](Particle &p1, Particle &p2) {
            double delta_x = p1.getX()[0] - p2.getX()[0];
            double delta_y = p1.getX()[1] - p2.getX()[1];
            double scalar =
                    p1.getM() * p2.getM() * std::pow(1.0 / std::sqrt(delta_x * delta_x + delta_y * delta_y), 3);
            double F_X = -delta_x * scalar;
            double F_Y = -delta_y * scalar;
            p1.add_to_F({F_X, F_Y, 0.});
            p2.add_to_F({-F_X, -F_Y, 0.});
        };
        fpairFun = fastPairFunction;
    }

    void FGravity::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
        setPairFun();
    }

    fpair_fun_t FGravity::getFastForceFunction() {
        return fpairFun;
    }

    static void fastForeAlt(std::vector<double> &force,
                     std::vector<double> &x,
                     std::vector<double> &eps,
                     std::vector<double> &sig,
                     std::vector<double> &m,
                     std::vector<int> &t,
                     unsigned long indexI,
                     double xJ0, double xJ1, double xJ2,
                     double, double, double mJ, int) {
        double d0, d1, d2, s, df0, df1, df2;
        d0 = x[3*indexI + 0] - xJ0;
        d1 = x[3*indexI + 1] - xJ1;
        d2 = x[3*indexI + 2] - xJ2;
        s = m[indexI] * mJ * std::pow(1.0 / std::sqrt(d0 * d0 + d1 * d1 + d2 * d2), 3);
        df0 = d0 * s;
        df1 = d1 * s;
        df2 = d2 * s;
        force[3*indexI + 0] -= df0;
        force[3*indexI + 1] -= df1;
        force[3*indexI + 2] -= df2;
    }

    fpair_fun_alt_t FGravity::getFastForceAltFunction() {
        return fastForeAlt;
    }

    static std::array<double,3> fastForceRet(std::vector<double> &force,
                             std::vector<double> &x,
                             std::vector<double> &eps,
                             std::vector<double> &sig,
                             std::vector<double> &m,
                             std::vector<int> &t,
                             unsigned long indexI,
                             double xJ0, double xJ1, double xJ2,
                             double, double, double mJ, int){
        double d0, d1, d2, s, df0, df1, df2;
        d0 = x[3*indexI + 0] - xJ0;
        d1 = x[3*indexI + 1] - xJ1;
        d2 = x[3*indexI + 2] - xJ2;
        s = m[indexI] * mJ * std::pow(1.0 / std::sqrt(d0 * d0 + d1 * d1 + d2 * d2), 3);
        df0 = d0 * s;
        df1 = d1 * s;
        df2 = d2 * s;
        return {-df0, -df1, -df2};
    }

    fpair_fun_ret_t FGravity::getFastForceRetFunction() {
        return fastForceRet;
    }
} // sim::physics::force