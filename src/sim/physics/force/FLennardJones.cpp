//
// Created by alex on 26.11.2022.
//

#include "FLennardJones.h"

namespace sim::physics::force {
    /**
     * @brief Calculates the using the All-Pairs algorithm. As displayed in the plot the Linked-Cell algorithm is much more performant in comparison
     * \image html plot.png width=800px
     * \image latex plot.eps "Runtime comparison of All-Pairs algorithm with Linked-Cell algorithm" width=5cm
     */
    void FLennardJones::operator()() {
        particleContainer.runOnActiveData([this](
                                       Kokkos::View<double*> &force,
                                       Kokkos::View<double*> &oldForce,
                                       Kokkos::View<double*> &x,
                                       Kokkos::View<double*> &v,
                                       Kokkos::View<double*> &m,
                                       Kokkos::View<int*> &type,
                                       unsigned long count,
                                       Kokkos::View<double*> &eps,
                                       Kokkos::View<double*> &sig,
                                       std::vector<unsigned long>& activeParticles){
            for(unsigned long indexI = 0; indexI < activeParticles.size(); indexI++){
                for(unsigned long indexJ = indexI + 1; indexJ < activeParticles.size(); indexJ++) {
                    this->fpairFun(force, x, eps, sig, m, type, activeParticles[indexI], activeParticles[indexJ]);
                }
            }
        });
    }

    pair_fun_t &FLennardJones::getForceFunction() {
        return pairFun;
    }

    void FLennardJones::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
        setPairFun();
    }

    static const double rt3_2 = std::pow(2,1.0/3.0);
    static void fastPairFunction(Kokkos::View<double*> &force,
                                 Kokkos::View<double*> &x,
                                 Kokkos::View<double*> &eps,
                                 Kokkos::View<double*> &sig,
                                 Kokkos::View<double*> &m,
                                 Kokkos::View<int*> &t,
                                 unsigned long indexI, unsigned long indexJ) {
        double sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1;
        sigma = (sig[indexI] + sig[indexJ]) / 2;
        sigma2 = sigma * sigma;
        sigma6 = sigma2 * sigma2 * sigma2;
        epsilon = std::sqrt(eps[indexI] * eps[indexJ]); // TODO this can be cached
        d0 = x[indexI*3 + 0] - x[indexJ*3 + 0];
        d1 = x[indexI*3 + 1] - x[indexJ*3 + 1];
        d2 = x[indexI*3 + 2] - x[indexJ*3 + 2];
        dsqr = d0*d0 + d1*d1 + d2*d2;
        //check if is membrane -> need to skip attractive forces
        if (t[indexI] & 0x80000000 || t[indexJ] & 0x80000000) {
            if (dsqr >= rt3_2 * sigma2) return;
        }

        l2NInvSquare = 1 / (dsqr);
        fac0 = 24 * epsilon * l2NInvSquare;
        l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;
        fac1_sum1 = sigma6 * l2NInvPow6;
        fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);

        force[indexI*3 + 0] -= fac0 * fac1 * d0;
        force[indexI*3 + 1] -= fac0 * fac1 * d1;
        force[indexI*3 + 2] -= fac0 * fac1 * d2;
        force[indexJ*3 + 0] += fac0 * fac1 * d0;
        force[indexJ*3 + 1] += fac0 * fac1 * d1;
        force[indexJ*3 + 2] += fac0 * fac1 * d2;
    }

    void FLennardJones::setPairFun() {
        pairFun = [](Particle &p1, Particle &p2) {
            Eigen::Vector3d delta{p1.getX() - p2.getX()};
            double l2Norm = delta.squaredNorm();
            double l2NInvSquare = 1 / (l2Norm * l2Norm);                        // invert squared norm
            double epsilon = std::sqrt(p1.getEpsilon()*p2.getEpsilon());
            double fac0 = 24 * epsilon * l2NInvSquare;                          // create first factor
            double l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;     // sixth power of inverted l2 norm
            double sigma = (p1.getSigma()+p2.getSigma())/2;
            double sigma6 = sigma * sigma * sigma;
            sigma6 = sigma6 * sigma6;                                           // sixth power of sigma
            double fac1_sum1 = sigma6 * l2NInvPow6;                             // first summand of middle factor
            double fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);            // create middle factor

            Eigen::Vector3d force{(-1) * fac0 * fac1 * delta};                  // bring it all together
            p1.add_to_F(force);
            p2.add_to_F(-force);                                                // reuse fact that F_ij = -F_ji
        };
        fpairFun = fastPairFunction;
    }

    fpair_fun_t FLennardJones::getFastForceFunction() {
        return fpairFun;
    }

    static void fastPairAlt(Kokkos::View<double*> &force,
                            Kokkos::View<double*> &x,
                            Kokkos::View<double*> &eps,
                            Kokkos::View<double*> &sig,
                            Kokkos::View<double*> &m,
                            Kokkos::View<int*> &t,
                            unsigned long indexI,
                            double xJ0, double xJ1, double xJ2,
                            double epsJ, double sigJ, double mJ, int tJ) {
        double sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1;
        sigma = (sig[indexI] + sigJ) / 2;
        sigma2 = sigma * sigma;
        sigma6 = sigma2 * sigma2 * sigma2;
        epsilon = std::sqrt(eps[indexI] * epsJ); // TODO this can be cached
        d0 = x[indexI*3 + 0] - xJ0;
        d1 = x[indexI*3 + 1] - xJ1;
        d2 = x[indexI*3 + 2] - xJ2;
        dsqr = d0*d0 + d1*d1 + d2*d2;
        //check if is membrane -> need to skip attractive forces
        if (t[indexI] & 0x80000000 || tJ & 0x80000000) {
            if (dsqr >= rt3_2 * sigma2) return;
        }

        l2NInvSquare = 1 / (dsqr);
        fac0 = 24 * epsilon * l2NInvSquare;
        l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;
        fac1_sum1 = sigma6 * l2NInvPow6;
        fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);

        force[indexI*3 + 0] -= fac0 * fac1 * d0;
        force[indexI*3 + 1] -= fac0 * fac1 * d1;
        force[indexI*3 + 2] -= fac0 * fac1 * d2;
    }

    fpair_fun_alt_t FLennardJones::getFastForceAltFunction() {
        return fastPairAlt;
    }

    std::array<double,3> fastPairRet(Kokkos::View<double*> &force,
                                     Kokkos::View<double*> &x,
                                     Kokkos::View<double*> &eps,
                                     Kokkos::View<double*> &sig,
                                     Kokkos::View<double*> &m,
                                     Kokkos::View<int*> &t,
                                     unsigned long indexI,
                                     double xJ0, double xJ1, double xJ2,
                                     double epsJ, double sigJ, double mJ, int tJ){
        double sigma, sigma2, sigma6, epsilon, d0, d1, d2, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1;
        sigma = (sig[indexI] + sigJ) / 2;
        sigma2 = sigma * sigma;
        sigma6 = sigma2 * sigma2 * sigma2;
        epsilon = std::sqrt(eps[indexI] * epsJ); // TODO this can be cached
        d0 = x[indexI*3 + 0] - xJ0;
        d1 = x[indexI*3 + 1] - xJ1;
        d2 = x[indexI*3 + 2] - xJ2;
        dsqr = d0*d0 + d1*d1 + d2*d2;
        //check if is membrane -> need to skip attractive forces
        if (t[indexI] & 0x80000000 || tJ & 0x80000000) {
            if (dsqr >= rt3_2 * sigma2) return {0,0,0};
        }

        l2NInvSquare = 1 / (dsqr);
        fac0 = 24 * epsilon * l2NInvSquare;
        l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;
        fac1_sum1 = sigma6 * l2NInvPow6;
        fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);

        return {fac0 * fac1 * d0,
                fac0 * fac1 * d1,
                fac0 * fac1 * d2};
    }

    fpair_fun_ret_t FLennardJones::getFastForceRetFunction() {
        return fastPairRet;
    }
} // sim::physics::force
