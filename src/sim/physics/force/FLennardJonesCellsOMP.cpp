//
// Created by alex on 11.01.23.
//

#include "FLennardJonesCellsOMP.h"
#include "defaults.h"

#include <iostream>
#include <immintrin.h>

namespace sim::physics::force {
    /**
     * @brief Returns the force function used
     *
     * @return pair_fun_t&
     */
    pair_fun_t &FLennardJonesCellsOMP::getForceFunction() {
        return pairFun;
    }

    /**
     * @brief The name says it all
     *
     * @param pc
     */
    void FLennardJonesCellsOMP::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
        forceDelegate.setParticleContainer(pc);
        setPairFun();
    }

    void FLennardJonesCellsOMP::setPairFun() {
        pairFun = forceDelegate.getForceFunction();
        fpairFun = forceDelegate.getFastForceFunction();
        fpairFunAlt = forceDelegate.getFastForceAltFunction();
        fpairFunRet = forceDelegate.getFastForceRetFunction();
    }

    fpair_fun_ret_t FLennardJonesCellsOMP::getFastForceRetFunction() {
        return fpairFunRet;
    }

    fpair_fun_t FLennardJonesCellsOMP::getFastForceFunction() {
        return fpairFun;
    }

    fpair_fun_alt_t FLennardJonesCellsOMP::getFastForceAltFunction() {
        return fpairFunAlt;
    }

    void FLennardJonesCellsOMP::operator()() {
        particleContainer.runOnDataCell([&](vec4d_t &force,
                                            vec4d_t &oldForce,
                                            vec4d_t &x,
                                            vec4d_t &v,
                                            std::vector<double> &m,
                                            std::vector<int> &type,
                                            unsigned long count,
                                            ParticleContainer::VectorCoordWrapper& cells,
                                            std::vector<double> &eps,
                                            std::vector<double> &sig){
            auto fpairFun = this->fpairFun;
            #pragma omp parallel for default(none) shared(cells, x, eps, sig, m, type, force, fpairFun) //reduction(+:interactions)
            for(size_t cellIndex=0; cellIndex < cells.size(); cellIndex++){
                auto& cell = cells[cellIndex];
                for(size_t i = 0; i < cell.size(); i++){
                    for(size_t j = i+1; j < cell.size(); j++){
                        fpairFun(force, x, eps, sig, m, type, cell[i], cell[j]);
                    }
                }
            }
        });

        particleContainer.runOnDataCell([&](vec4d_t &force,
                                            vec4d_t &oldForce,
                                            vec4d_t &x,
                                            vec4d_t &v,
                                            std::vector<double> &m,
                                            std::vector<int> &t,
                                            unsigned long count,
                                            ParticleContainer::VectorCoordWrapper &cells,
                                            std::vector<double> &eps,
                                            std::vector<double> &sig) {
            static const double rt3_2 = std::pow(2, 1 / 3);
            const std::vector<std::vector<std::pair<unsigned long, unsigned long>>>& alternativeTaskGroups = particleContainer.generateDistinctAlternativeCellNeighbours();
            double *_force = force.data();
            double *_x = x.data();
            double *_sig = sig.data();
            double *_eps = eps.data();
            size_t size = force.size();
            double sigma, sigma2, sigma6, epsilon, dsqr, l2NInvSquare, fac0, l2NInvPow6, fac1_sum1, fac1;
            unsigned long indexI;
            unsigned long indexJ;
            unsigned long indexII;
            unsigned long indexJJ;
            unsigned long indexX;
            unsigned long indexY;
            unsigned long indexC0;
            unsigned long indexC1;
            __m256d* _vf = reinterpret_cast<__m256d*>(_force);
            size_t vf_size = size/4;

            size_t maxThreads = omp_get_max_threads();
            #pragma omp declare reduction \
                (addpd:__m256d:omp_out+=omp_in) initializer(omp_priv=_mm256_setzero_pd())

            #pragma omp parallel \
                default(none) \
                shared(size, _x, x, t, cells, _eps, eps, _sig, sig,alternativeTaskGroups,fpairFun,force,m, maxThreads, vf_size) \
                private(sigma, sigma2, sigma6, epsilon, dsqr, l2NInvSquare, \
                        fac0, l2NInvPow6, fac1_sum1, fac1, indexI, indexJ,indexII,indexJJ,indexY,indexC0,indexC1) \
                firstprivate(indexX, rt3_2) \
                reduction(addpd:_vf[:vf_size])
            {
                //generate tasks: for all distinct cell neighbours
                #pragma omp for
                for(indexX = 0; indexX < maxThreads; indexX++) {
                    indexX = omp_get_thread_num();
                    for(indexY = 0; indexY < alternativeTaskGroups[indexX].size(); indexY++){
                        indexC0 = alternativeTaskGroups[indexX][indexY].first;
                        indexC1 = alternativeTaskGroups[indexX][indexY].second;
                        //TASK:
                        //first do 0-3 iterations for cellI s.t. its size%4==0
                        //handle cellJ: do make size%4==0, then vectorize
                        //then handle rest of cellI vectorized
                        //handle cellJ in the same way

                        //intentionally doing loop unrolling -> more performance but more code
                        //init
                        indexII = 0;
                        static const __m256i xMask = _mm256_set_epi64x(0, -1, -1, -1);
                        static const __m256d half = _mm256_set1_pd(0.5);
                        static const __m256d one = _mm256_set1_pd(1.0);
                        static const __m256d two = _mm256_set1_pd(2.0);
                        static const __m256d fac24 = _mm256_set1_pd(24.0);

                        // make cellI size divisible by 4
                        for(; indexII < cells[indexC0].size() % 4; indexII++) {
                            indexJJ = 0;
                            indexI = cells[indexC0][indexII];
                            double sigIs = sig[indexI];
                            double epsIs = eps[indexI];
                            __m256d xI = _mm256_maskload_pd(_x + 4 * indexI, xMask);
                            __m256d sigI = _mm256_set1_pd(sigIs); //sigma of I in all positions
                            __m256d epsI = _mm256_set1_pd(epsIs); //epsilon of I in all positions

                            //I is not vector, J is not vector
                            for(; indexJJ < cells[indexC1].size() % 4; indexJJ++) {
                                indexJ = cells[indexC1][indexJJ];
                                sigma = (sigIs + sig[indexJ]) / 2;
                                epsilon = std::sqrt(epsIs * eps[indexJ]);
                                sigma6 = std::pow(sigma,6);

                                __m256d xJ = _mm256_maskload_pd(_x + 4 * indexJ, xMask);
                                __m256d d = _mm256_sub_pd(xI, xJ);
                                __m256d d_sqr = _mm256_mul_pd(d,d);
                                __m256d hadd = _mm256_hadd_pd(d_sqr, d_sqr); //dq0+dq1 -> [63:0], dq2+0 -> [191:128]
                                __m128d upper = _mm256_extractf128_pd(hadd, 1);
                                __m128d d_sum = _mm_add_pd(upper, _mm256_castpd256_pd128(hadd));
                                dsqr = _mm_cvtsd_f64(d_sum);

                                l2NInvSquare = 1 / (dsqr);
                                fac0 = 24 * epsilon * l2NInvSquare;
                                l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;
                                fac1_sum1 = sigma6 * l2NInvPow6;
                                fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);
                                __m256d scalar = _mm256_set1_pd(fac0 * fac1);
                                __m256d result = _mm256_mul_pd(scalar, d);
                                _vf[indexI] -= result;
                                _vf[indexJ] += result;
                            }
                            //I is not vector, J is vector now
                            for(; indexJJ < cells[indexC1].size(); indexJJ += 4) {
                                indexJ = cells[indexC1][indexJJ];
                                __m256d sigJ = _mm256_loadu_pd(_sig + indexJ); // sigma of all 4 particles
                                __m256d epsJ = _mm256_loadu_pd(_eps + indexJ); // epsilon of all 4 particles
                                __m256d xJ0 = _mm256_maskload_pd(_x + 4 * indexJ, xMask);
                                __m256d xJ1 = _mm256_maskload_pd(_x + 4 * indexJ + 4, xMask);
                                __m256d xJ2 = _mm256_maskload_pd(_x + 4 * indexJ + 8, xMask);
                                __m256d xJ3 = _mm256_maskload_pd(_x + 4 * indexJ + 12, xMask);

                                __m256d tmpSig = _mm256_add_pd(sigI, sigJ);
                                __m256d sig = _mm256_mul_pd(tmpSig, half);
                                __m256d tmpEps = _mm256_mul_pd(epsI, epsJ);
                                __m256d eps = _mm256_sqrt_pd(tmpEps);

                                __m256d d0 = _mm256_sub_pd(xI, xJ0);
                                __m256d d1 = _mm256_sub_pd(xI, xJ1);
                                __m256d d2 = _mm256_sub_pd(xI, xJ2);
                                __m256d d3 = _mm256_sub_pd(xI, xJ3);
                                __m256d d0_sqr = _mm256_mul_pd(d0,d0);
                                __m256d d1_sqr = _mm256_mul_pd(d1,d1);
                                __m256d d2_sqr = _mm256_mul_pd(d2,d2);
                                __m256d d3_sqr = _mm256_mul_pd(d3,d3);
                                __m256d hadd10 = _mm256_hadd_pd(d0_sqr, d1_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                __m256d hadd32 = _mm256_hadd_pd(d2_sqr, d3_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                __m128d upper10 = _mm256_extractf128_pd(hadd10, 1);
                                __m128d upper32 = _mm256_extractf128_pd(hadd32, 1);
                                __m128d d10 = _mm_add_pd(upper10, _mm256_castpd256_pd128(hadd10));
                                __m128d d32 = _mm_add_pd(upper32, _mm256_castpd256_pd128(hadd32));
                                __m256d d3210 = _mm256_set_m128d(d32, d10); // this is dsqr for all 4 particles

                                __m256d l2inv_dsqr_3210 = _mm256_div_pd(one, d3210);
                                __m256d fac0_3210 = _mm256_mul_pd(fac24, eps);
                                        fac0_3210 = _mm256_mul_pd(fac0_3210, l2inv_dsqr_3210);
                                __m256d l2inv_pow6_3210 = _mm256_mul_pd(l2inv_dsqr_3210, _mm256_mul_pd(l2inv_dsqr_3210,l2inv_dsqr_3210));
                                __m256d sigP2 = _mm256_mul_pd(sig, sig);
                                __m256d sigP6 = _mm256_mul_pd(sigP2, _mm256_mul_pd(sigP2, sigP2));
                                __m256d fac1_tmp_3210 = _mm256_mul_pd(sigP6, l2inv_pow6_3210);
                                __m256d fac1_3210 = _mm256_sub_pd(fac1_tmp_3210, _mm256_mul_pd(two, _mm256_mul_pd(fac1_tmp_3210,fac1_tmp_3210)));
                                __m256d scale_3210 = _mm256_mul_pd(fac0_3210, fac1_3210);
                                __m128d scale_32 = _mm256_extractf128_pd(scale_3210, 1);
                                __m128d scale_10  = _mm256_castpd256_pd128(scale_3210);
                                __m256d scale_0 = _mm256_broadcastsd_pd(scale_10);
                                __m256d scale_1 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_10, 1));
                                __m256d scale_2 = _mm256_broadcastsd_pd(scale_32);
                                __m256d scale_3 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_32, 1));
                                __m256d res0 = _mm256_mul_pd(scale_0, d0);
                                __m256d res1 = _mm256_mul_pd(scale_1, d1);
                                __m256d res2 = _mm256_mul_pd(scale_2, d2);
                                __m256d res3 = _mm256_mul_pd(scale_3, d3);
                                _vf[indexI] -= res0;
                                _vf[indexI] -= res1;
                                _vf[indexI] -= res2;
                                _vf[indexI] -= res3;
                                _vf[indexJ + 0] += res0;
                                _vf[indexJ + 1] += res1;
                                _vf[indexJ + 2] += res2;
                                _vf[indexJ + 3] += res3;
                            }
                        }
                        // cellI is divisible by 4 now
                        for(; indexII < cells[indexC0].size(); indexII += 4){
                            indexJJ = 0;
                            indexI = cells[indexC0][indexII];
                            __m256d sigI = _mm256_loadu_pd(_sig + indexI); // 4 diff. sigma values for I
                            __m256d epsI = _mm256_loadu_pd(_eps + indexI); // 4 diff. epsilon values for I
                            __m256d xI0 = _mm256_maskload_pd(_x + 4 * indexI, xMask);
                            __m256d xI1 = _mm256_maskload_pd(_x + 4 * indexI + 4, xMask);
                            __m256d xI2 = _mm256_maskload_pd(_x + 4 * indexI + 8, xMask);
                            __m256d xI3 = _mm256_maskload_pd(_x + 4 * indexI + 12, xMask);

                            //I is vector, J is not necessarily vector
                            for(; indexJJ < cells[indexC1].size() % 4; indexJJ++) {
                                indexJ = cells[indexC1][indexJJ];
                                __m256d sigJ = _mm256_set1_pd(_sig[indexJ]);
                                __m256d epsJ = _mm256_set1_pd(_eps[indexJ]);
                                __m256d xJ = _mm256_maskload_pd(_x + 4 * indexJ, xMask);

                                __m256d tmpSig = _mm256_add_pd(sigI, sigJ);
                                __m256d sig = _mm256_mul_pd(tmpSig, half);
                                __m256d tmpEps = _mm256_mul_pd(epsI, epsJ);
                                __m256d eps = _mm256_sqrt_pd(tmpEps);

                                __m256d d0 = _mm256_sub_pd(xI0, xJ);
                                __m256d d1 = _mm256_sub_pd(xI1, xJ);
                                __m256d d2 = _mm256_sub_pd(xI2, xJ);
                                __m256d d3 = _mm256_sub_pd(xI3, xJ);
                                __m256d d0_sqr = _mm256_mul_pd(d0,d0);
                                __m256d d1_sqr = _mm256_mul_pd(d1,d1);
                                __m256d d2_sqr = _mm256_mul_pd(d2,d2);
                                __m256d d3_sqr = _mm256_mul_pd(d3,d3);
                                __m256d hadd10 = _mm256_hadd_pd(d0_sqr, d1_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                __m256d hadd32 = _mm256_hadd_pd(d2_sqr, d3_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                __m128d upper10 = _mm256_extractf128_pd(hadd10, 1);
                                __m128d upper32 = _mm256_extractf128_pd(hadd32, 1);
                                __m128d d10 = _mm_add_pd(upper10, _mm256_castpd256_pd128(hadd10));
                                __m128d d32 = _mm_add_pd(upper32, _mm256_castpd256_pd128(hadd32));
                                __m256d d3210 = _mm256_set_m128d(d32, d10); // this is dsqr for all 4 particles

                                __m256d l2inv_dsqr_3210 = _mm256_div_pd(one, d3210);
                                __m256d fac0_3210 = _mm256_mul_pd(fac24, eps);
                                fac0_3210 = _mm256_mul_pd(fac0_3210, l2inv_dsqr_3210);
                                __m256d l2inv_pow6_3210 = _mm256_mul_pd(l2inv_dsqr_3210, _mm256_mul_pd(l2inv_dsqr_3210,l2inv_dsqr_3210));
                                __m256d sigP2 = _mm256_mul_pd(sig, sig);
                                __m256d sigP6 = _mm256_mul_pd(sigP2, _mm256_mul_pd(sigP2, sigP2));
                                __m256d fac1_tmp_3210 = _mm256_mul_pd(sigP6, l2inv_pow6_3210);
                                __m256d fac1_3210 = _mm256_sub_pd(fac1_tmp_3210, _mm256_mul_pd(two, _mm256_mul_pd(fac1_tmp_3210,fac1_tmp_3210)));
                                __m256d scale_3210 = _mm256_mul_pd(fac0_3210, fac1_3210);
                                __m128d scale_32 = _mm256_extractf128_pd(scale_3210, 1);
                                __m128d scale_10  = _mm256_castpd256_pd128(scale_3210);
                                __m256d scale_0 = _mm256_broadcastsd_pd(scale_10);
                                __m256d scale_1 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_10, 1));
                                __m256d scale_2 = _mm256_broadcastsd_pd(scale_32);
                                __m256d scale_3 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_32, 1));
                                __m256d res0 = _mm256_mul_pd(scale_0, d0);
                                __m256d res1 = _mm256_mul_pd(scale_1, d1);
                                __m256d res2 = _mm256_mul_pd(scale_2, d2);
                                __m256d res3 = _mm256_mul_pd(scale_3, d3);
                                _vf[indexI + 0] -= res0;
                                _vf[indexI + 1] -= res1;
                                _vf[indexI + 1] -= res2;
                                _vf[indexI + 3] -= res3;
                                _vf[indexJ] += res0;
                                _vf[indexJ] += res1;
                                _vf[indexJ] += res2;
                                _vf[indexJ] += res3;
                            }
                            //I is vector, J is vector -> need to let 4x4 particles interact with each other
                            for(; indexJJ < cells[indexC1].size(); indexJJ += 4) {
                                indexJ = cells[indexC1][indexJJ];
                                __m256d sigJ = _mm256_loadu_pd(_sig + indexJ);
                                __m256d epsJ = _mm256_loadu_pd(_eps + indexJ);
                                //rotate sigI right etc and generate 16 diff. pairs

                                __m256d xJ0 = _mm256_maskload_pd(_x + 4 * indexJ, xMask);
                                __m256d xJ1 = _mm256_maskload_pd(_x + 4 * indexJ + 4, xMask);
                                __m256d xJ2 = _mm256_maskload_pd(_x + 4 * indexJ + 8, xMask);
                                __m256d xJ3 = _mm256_maskload_pd(_x + 4 * indexJ + 12, xMask);

                                //handle different pairs sequentially in blocks of 4 to not overuse registers
                                //I0
                                {
                                    __m256d dI0J0 = _mm256_sub_pd(xI0, xJ0);
                                    __m256d dI0J1 = _mm256_sub_pd(xI0, xJ1);
                                    __m256d dI0J2 = _mm256_sub_pd(xI0, xJ2);
                                    __m256d dI0J3 = _mm256_sub_pd(xI0, xJ3);
                                    __m256d sigI_tmp = _mm256_set1_pd(_mm256_cvtsd_f64(sigI));
                                    __m256d sig_tmp = _mm256_add_pd(sigI_tmp, sigJ);
                                    __m256d sig = _mm256_mul_pd(sig_tmp, half);
                                    __m256d epsI_tmp = _mm256_set1_pd(_mm256_cvtsd_f64(epsI));
                                    __m256d eps_tmp = _mm256_mul_pd(epsI_tmp, epsJ);
                                    __m256d eps = _mm256_sqrt_pd(eps_tmp);

                                    __m256d dI0J0_sqr = _mm256_mul_pd(dI0J0, dI0J0); //0
                                    __m256d dI0J1_sqr = _mm256_mul_pd(dI0J1, dI0J1); //1
                                    __m256d dI0J2_sqr = _mm256_mul_pd(dI0J2, dI0J2); //2
                                    __m256d dI0J3_sqr = _mm256_mul_pd(dI0J3, dI0J3); //3
                                    __m256d hadd10 = _mm256_hadd_pd(dI0J0_sqr, dI0J1_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                    __m256d hadd32 = _mm256_hadd_pd(dI0J2_sqr, dI0J3_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                    __m128d upper10 = _mm256_extractf128_pd(hadd10, 1);
                                    __m128d upper32 = _mm256_extractf128_pd(hadd32, 1);
                                    __m128d d10 = _mm_add_pd(upper10, _mm256_castpd256_pd128(hadd10));
                                    __m128d d32 = _mm_add_pd(upper32, _mm256_castpd256_pd128(hadd32));
                                    __m256d d3210 = _mm256_set_m128d(d32, d10); // this is dsqr for all 4 particles

                                    __m256d l2inv_dsqr_3210 = _mm256_div_pd(one, d3210);
                                    __m256d fac0_3210 = _mm256_mul_pd(fac24, eps);
                                    fac0_3210 = _mm256_mul_pd(fac0_3210, l2inv_dsqr_3210);
                                    __m256d l2inv_pow6_3210 = _mm256_mul_pd(l2inv_dsqr_3210, _mm256_mul_pd(l2inv_dsqr_3210,l2inv_dsqr_3210));
                                    __m256d sigP2 = _mm256_mul_pd(sig, sig);
                                    __m256d sigP6 = _mm256_mul_pd(sigP2, _mm256_mul_pd(sigP2, sigP2));
                                    __m256d fac1_tmp_3210 = _mm256_mul_pd(sigP6, l2inv_pow6_3210);
                                    __m256d fac1_3210 = _mm256_sub_pd(fac1_tmp_3210, _mm256_mul_pd(two, _mm256_mul_pd(fac1_tmp_3210,fac1_tmp_3210)));
                                    __m256d scale_3210 = _mm256_mul_pd(fac0_3210, fac1_3210);
                                    __m128d scale_32 = _mm256_extractf128_pd(scale_3210, 1);
                                    __m128d scale_10  = _mm256_castpd256_pd128(scale_3210);
                                    __m256d scale_0 = _mm256_broadcastsd_pd(scale_10);
                                    __m256d scale_1 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_10, 1));
                                    __m256d scale_2 = _mm256_broadcastsd_pd(scale_32);
                                    __m256d scale_3 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_32, 1));
                                    __m256d res0 = _mm256_mul_pd(scale_0, dI0J0);
                                    __m256d res1 = _mm256_mul_pd(scale_1, dI0J1);
                                    __m256d res2 = _mm256_mul_pd(scale_2, dI0J2);
                                    __m256d res3 = _mm256_mul_pd(scale_3, dI0J3);
                                    _vf[indexI] -= res0;
                                    _vf[indexI] -= res1;
                                    _vf[indexI] -= res2;
                                    _vf[indexI] -= res3;
                                    _vf[indexJ + 0] += res0;
                                    _vf[indexJ + 1] += res1;
                                    _vf[indexJ + 2] += res2;
                                    _vf[indexJ + 3] += res3;
                                }

                                //I1
                                {
                                    __m256d dI1J0 = _mm256_sub_pd(xI1, xJ0);
                                    __m256d dI1J1 = _mm256_sub_pd(xI1, xJ1);
                                    __m256d dI1J2 = _mm256_sub_pd(xI1, xJ2);
                                    __m256d dI1J3 = _mm256_sub_pd(xI1, xJ3);
                                    __m256d sigI_tmp = _mm256_set1_pd(_mm256_cvtsd_f64(_mm256_permute_pd(sigI, 1)));
                                    __m256d sig_tmp = _mm256_add_pd(sigI_tmp, sigJ);
                                    __m256d sig = _mm256_mul_pd(sig_tmp, half);
                                    __m256d epsI_tmp = _mm256_set1_pd(_mm256_cvtsd_f64(_mm256_permute_pd(epsI, 1)));
                                    __m256d eps_tmp = _mm256_mul_pd(epsI_tmp, epsJ);
                                    __m256d eps = _mm256_sqrt_pd(eps_tmp);

                                    __m256d dI1J0_sqr = _mm256_mul_pd(dI1J0, dI1J0); //0
                                    __m256d dI1J1_sqr = _mm256_mul_pd(dI1J1, dI1J1); //1
                                    __m256d dI1J2_sqr = _mm256_mul_pd(dI1J2, dI1J2); //2
                                    __m256d dI1J3_sqr = _mm256_mul_pd(dI1J3, dI1J3); //3
                                    __m256d hadd10 = _mm256_hadd_pd(dI1J0_sqr, dI1J1_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                    __m256d hadd32 = _mm256_hadd_pd(dI1J2_sqr, dI1J3_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                    __m128d upper10 = _mm256_extractf128_pd(hadd10, 1);
                                    __m128d upper32 = _mm256_extractf128_pd(hadd32, 1);
                                    __m128d d10 = _mm_add_pd(upper10, _mm256_castpd256_pd128(hadd10));
                                    __m128d d32 = _mm_add_pd(upper32, _mm256_castpd256_pd128(hadd32));
                                    __m256d d3210 = _mm256_set_m128d(d32, d10); // this is dsqr for all 4 particles

                                    __m256d l2inv_dsqr_3210 = _mm256_div_pd(one, d3210);
                                    __m256d fac0_3210 = _mm256_mul_pd(fac24, eps);
                                    fac0_3210 = _mm256_mul_pd(fac0_3210, l2inv_dsqr_3210);
                                    __m256d l2inv_pow6_3210 = _mm256_mul_pd(l2inv_dsqr_3210, _mm256_mul_pd(l2inv_dsqr_3210,l2inv_dsqr_3210));
                                    __m256d sigP2 = _mm256_mul_pd(sig, sig);
                                    __m256d sigP6 = _mm256_mul_pd(sigP2, _mm256_mul_pd(sigP2, sigP2));
                                    __m256d fac1_tmp_3210 = _mm256_mul_pd(sigP6, l2inv_pow6_3210);
                                    __m256d fac1_3210 = _mm256_sub_pd(fac1_tmp_3210, _mm256_mul_pd(two, _mm256_mul_pd(fac1_tmp_3210,fac1_tmp_3210)));
                                    __m256d scale_3210 = _mm256_mul_pd(fac0_3210, fac1_3210);
                                    __m128d scale_32 = _mm256_extractf128_pd(scale_3210, 1);
                                    __m128d scale_10  = _mm256_castpd256_pd128(scale_3210);
                                    __m256d scale_0 = _mm256_broadcastsd_pd(scale_10);
                                    __m256d scale_1 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_10, 1));
                                    __m256d scale_2 = _mm256_broadcastsd_pd(scale_32);
                                    __m256d scale_3 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_32, 1));
                                    __m256d res0 = _mm256_mul_pd(scale_0, dI1J0);
                                    __m256d res1 = _mm256_mul_pd(scale_1, dI1J1);
                                    __m256d res2 = _mm256_mul_pd(scale_2, dI1J2);
                                    __m256d res3 = _mm256_mul_pd(scale_3, dI1J3);
                                    _vf[indexI + 1] -= res0;
                                    _vf[indexI + 1] -= res1;
                                    _vf[indexI + 1] -= res2;
                                    _vf[indexI + 1] -= res3;
                                    _vf[indexJ + 0] += res0;
                                    _vf[indexJ + 1] += res1;
                                    _vf[indexJ + 2] += res2;
                                    _vf[indexJ + 3] += res3;
                                }

                                //I2
                                {
                                    __m256d dI2J0 = _mm256_sub_pd(xI2, xJ0);
                                    __m256d dI2J1 = _mm256_sub_pd(xI2, xJ1);
                                    __m256d dI2J2 = _mm256_sub_pd(xI2, xJ2);
                                    __m256d dI2J3 = _mm256_sub_pd(xI2, xJ3);
                                    __m256d sigI_tmp = _mm256_set1_pd(_mm_cvtsd_f64(_mm256_extractf128_pd(sigI, 1)));
                                    __m256d sig_tmp = _mm256_add_pd(sigI_tmp, sigJ);
                                    __m256d sig = _mm256_mul_pd(sig_tmp, half);
                                    __m256d epsI_tmp = _mm256_set1_pd(_mm_cvtsd_f64(_mm256_extractf128_pd(epsI, 1)));
                                    __m256d eps_tmp = _mm256_mul_pd(epsI_tmp, epsJ);
                                    __m256d eps = _mm256_sqrt_pd(eps_tmp);

                                    __m256d dI2J0_sqr = _mm256_mul_pd(dI2J0, dI2J0); //0
                                    __m256d dI2J1_sqr = _mm256_mul_pd(dI2J1, dI2J1); //1
                                    __m256d dI2J2_sqr = _mm256_mul_pd(dI2J2, dI2J2); //2
                                    __m256d dI2J3_sqr = _mm256_mul_pd(dI2J3, dI2J3); //3
                                    __m256d hadd10 = _mm256_hadd_pd(dI2J0_sqr, dI2J1_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                    __m256d hadd32 = _mm256_hadd_pd(dI2J2_sqr, dI2J3_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                    __m128d upper10 = _mm256_extractf128_pd(hadd10, 1);
                                    __m128d upper32 = _mm256_extractf128_pd(hadd32, 1);
                                    __m128d d10 = _mm_add_pd(upper10, _mm256_castpd256_pd128(hadd10));
                                    __m128d d32 = _mm_add_pd(upper32, _mm256_castpd256_pd128(hadd32));
                                    __m256d d3210 = _mm256_set_m128d(d32, d10); // this is dsqr for all 4 particles

                                    __m256d l2inv_dsqr_3210 = _mm256_div_pd(one, d3210);
                                    __m256d fac0_3210 = _mm256_mul_pd(fac24, eps);
                                    fac0_3210 = _mm256_mul_pd(fac0_3210, l2inv_dsqr_3210);
                                    __m256d l2inv_pow6_3210 = _mm256_mul_pd(l2inv_dsqr_3210, _mm256_mul_pd(l2inv_dsqr_3210,l2inv_dsqr_3210));
                                    __m256d sigP2 = _mm256_mul_pd(sig, sig);
                                    __m256d sigP6 = _mm256_mul_pd(sigP2, _mm256_mul_pd(sigP2, sigP2));
                                    __m256d fac1_tmp_3210 = _mm256_mul_pd(sigP6, l2inv_pow6_3210);
                                    __m256d fac1_3210 = _mm256_sub_pd(fac1_tmp_3210, _mm256_mul_pd(two, _mm256_mul_pd(fac1_tmp_3210,fac1_tmp_3210)));
                                    __m256d scale_3210 = _mm256_mul_pd(fac0_3210, fac1_3210);
                                    __m128d scale_32 = _mm256_extractf128_pd(scale_3210, 1);
                                    __m128d scale_10  = _mm256_castpd256_pd128(scale_3210);
                                    __m256d scale_0 = _mm256_broadcastsd_pd(scale_10);
                                    __m256d scale_1 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_10, 1));
                                    __m256d scale_2 = _mm256_broadcastsd_pd(scale_32);
                                    __m256d scale_3 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_32, 1));
                                    __m256d res0 = _mm256_mul_pd(scale_0, dI2J0);
                                    __m256d res1 = _mm256_mul_pd(scale_1, dI2J1);
                                    __m256d res2 = _mm256_mul_pd(scale_2, dI2J2);
                                    __m256d res3 = _mm256_mul_pd(scale_3, dI2J3);
                                    _vf[indexI + 2] -= res0;
                                    _vf[indexI + 2] -= res1;
                                    _vf[indexI + 2] -= res2;
                                    _vf[indexI + 2] -= res3;
                                    _vf[indexJ + 0] += res0;
                                    _vf[indexJ + 1] += res1;
                                    _vf[indexJ + 2] += res2;
                                    _vf[indexJ + 3] += res3;
                                }

                                //I3
                                {
                                    __m256d dI3J0 = _mm256_sub_pd(xI3, xJ0);
                                    __m256d dI3J1 = _mm256_sub_pd(xI3, xJ1);
                                    __m256d dI3J2 = _mm256_sub_pd(xI3, xJ2);
                                    __m256d dI3J3 = _mm256_sub_pd(xI3, xJ3);
                                    __m256d sigI_tmp = _mm256_set1_pd(_mm_cvtsd_f64(_mm_permute_pd(_mm256_extractf128_pd(sigI, 1),1)));
                                    __m256d sig_tmp = _mm256_add_pd(sigI_tmp, sigJ);
                                    __m256d sig = _mm256_mul_pd(sig_tmp, half);
                                    __m256d epsI_tmp = _mm256_set1_pd(_mm_cvtsd_f64(_mm_permute_pd(_mm256_extractf128_pd(epsI, 1),1)));
                                    __m256d eps_tmp = _mm256_mul_pd(epsI_tmp, epsJ);
                                    __m256d eps = _mm256_sqrt_pd(eps_tmp);

                                    __m256d dI3J0_sqr = _mm256_mul_pd(dI3J0, dI3J0); //0
                                    __m256d dI3J1_sqr = _mm256_mul_pd(dI3J1, dI3J1); //1
                                    __m256d dI3J2_sqr = _mm256_mul_pd(dI3J2, dI3J2); //2
                                    __m256d dI3J3_sqr = _mm256_mul_pd(dI3J3, dI3J3); //3
                                    __m256d hadd10 = _mm256_hadd_pd(dI3J0_sqr, dI3J1_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                    __m256d hadd32 = _mm256_hadd_pd(dI3J2_sqr, dI3J3_sqr); //(a,b) -> b3+b2, a3+a2, b1+b0, a1+a0
                                    __m128d upper10 = _mm256_extractf128_pd(hadd10, 1);
                                    __m128d upper32 = _mm256_extractf128_pd(hadd32, 1);
                                    __m128d d10 = _mm_add_pd(upper10, _mm256_castpd256_pd128(hadd10));
                                    __m128d d32 = _mm_add_pd(upper32, _mm256_castpd256_pd128(hadd32));
                                    __m256d d3210 = _mm256_set_m128d(d32, d10); // this is dsqr for all 4 particles

                                    __m256d l2inv_dsqr_3210 = _mm256_div_pd(one, d3210);
                                    __m256d fac0_3210 = _mm256_mul_pd(fac24, eps);
                                    fac0_3210 = _mm256_mul_pd(fac0_3210, l2inv_dsqr_3210);
                                    __m256d l2inv_pow6_3210 = _mm256_mul_pd(l2inv_dsqr_3210, _mm256_mul_pd(l2inv_dsqr_3210,l2inv_dsqr_3210));
                                    __m256d sigP2 = _mm256_mul_pd(sig, sig);
                                    __m256d sigP6 = _mm256_mul_pd(sigP2, _mm256_mul_pd(sigP2, sigP2));
                                    __m256d fac1_tmp_3210 = _mm256_mul_pd(sigP6, l2inv_pow6_3210);
                                    __m256d fac1_3210 = _mm256_sub_pd(fac1_tmp_3210, _mm256_mul_pd(two, _mm256_mul_pd(fac1_tmp_3210,fac1_tmp_3210)));
                                    __m256d scale_3210 = _mm256_mul_pd(fac0_3210, fac1_3210);
                                    __m128d scale_32 = _mm256_extractf128_pd(scale_3210, 1);
                                    __m128d scale_10  = _mm256_castpd256_pd128(scale_3210);
                                    __m256d scale_0 = _mm256_broadcastsd_pd(scale_10);
                                    __m256d scale_1 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_10, 1));
                                    __m256d scale_2 = _mm256_broadcastsd_pd(scale_32);
                                    __m256d scale_3 = _mm256_broadcastsd_pd(_mm_permute_pd(scale_32, 1));
                                    __m256d res0 = _mm256_mul_pd(scale_0, dI3J0);
                                    __m256d res1 = _mm256_mul_pd(scale_1, dI3J1);
                                    __m256d res2 = _mm256_mul_pd(scale_2, dI3J2);
                                    __m256d res3 = _mm256_mul_pd(scale_3, dI3J3);
                                    _vf[indexI + 3] -= res0;
                                    _vf[indexI + 3] -= res1;
                                    _vf[indexI + 3] -= res2;
                                    _vf[indexI + 3] -= res3;
                                    _vf[indexJ + 0] += res0;
                                    _vf[indexJ + 1] += res1;
                                    _vf[indexJ + 2] += res2;
                                    _vf[indexJ + 3] += res3;
                                }
                            }
                        }

//                        for (indexII = 0; indexII < cells[indexC0].size(); indexII++) {
//                            for (indexJJ = 0; indexJJ < cells[indexC1].size(); indexJJ++) {
//                                indexI = cells[indexC0][indexII];
//                                indexJ = cells[indexC1][indexJJ];
//                                //fpairFun(force, x, eps, sig, m, t, indexI, indexJ);
//                                sigma = (sig[indexI] + sig[indexJ]) / 2;
//                                sigma2 = sigma * sigma;
//                                sigma6 = sigma2 * sigma2 * sigma2;
//                                epsilon = std::sqrt(eps[indexI] * eps[indexJ]); // TODO this can be cached
//                                d0 = x[indexI * 3 + 0] - x[indexJ * 3 + 0];
//                                d1 = x[indexI * 3 + 1] - x[indexJ * 3 + 1];
//                                d2 = x[indexI * 3 + 2] - x[indexJ * 3 + 2];
//                                dsqr = d0 * d0 + d1 * d1 + d2 * d2;
//                                //check if is membrane -> need to skip attractive forces
//                                if (t[indexI] & 0x80000000 || t[indexJ] & 0x80000000) {
//                                    if (dsqr >= rt3_2 * sigma2) continue;
//                                }
//
//                                l2NInvSquare = 1 / (dsqr);
//                                fac0 = 24 * epsilon * l2NInvSquare;
//                                l2NInvPow6 = l2NInvSquare * l2NInvSquare * l2NInvSquare;
//                                fac1_sum1 = sigma6 * l2NInvPow6;
//                                fac1 = (fac1_sum1) - 2 * (fac1_sum1 * fac1_sum1);
//
//                                _force[indexI * 3 + 0] -= fac0 * fac1 * d0;
//                                _force[indexI * 3 + 1] -= fac0 * fac1 * d1;
//                                _force[indexI * 3 + 2] -= fac0 * fac1 * d2;
//                                _force[indexJ * 3 + 0] += fac0 * fac1 * d0;
//                                _force[indexJ * 3 + 1] += fac0 * fac1 * d1;
//                                _force[indexJ * 3 + 2] += fac0 * fac1 * d2;
//                            }
//                        }
                    }

                }
            }
        });
    }
} // force