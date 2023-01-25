//
// Created by alex on 11.01.23.
//

#include "FGravityCells.h"

namespace sim::physics::force {
    void FGravityCells::setPairFun() {
        pairFun = forceDelegate.getForceFunction();
        fpairFun = forceDelegate.getFastForceFunction();
        fpairFunAlt = forceDelegate.getFastForceAltFunction();
        fpairFunRet = forceDelegate.getFastForceRetFunction();
    }

    void FGravityCells::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
        forceDelegate.setParticleContainer(pc);
        setPairFun();
    }

    pair_fun_t &FGravityCells::getForceFunction() {
        return pairFun;
    }

    fpair_fun_t FGravityCells::getFastForceFunction() {
        return fpairFun;
    }

    void FGravityCells::operator()() {
        particleContainer.forAllCells([this](std::vector<double> &force,
                                             std::vector<double> &oldForce,
                                             std::vector<double> &x,
                                             std::vector<double> &v,
                                             std::vector<double> &m,
                                             std::vector<int> &type,
                                             unsigned long count,
                                             std::vector<unsigned long> &cellItems,
                                             std::vector<double> &eps,
                                             std::vector<double> &sig){
            for(unsigned long indexX = 0; indexX < cellItems.size(); indexX++){
                for(unsigned long indexY = indexX + 1; indexY < cellItems.size(); indexY++) {
                    unsigned long indexI = cellItems[indexX];
                    unsigned long indexJ = cellItems[indexY];
                    this->fpairFun(force, x, eps, sig, m, type, indexI, indexJ);
                }
            }
        });

        particleContainer.forAllDistinctCellNeighbours([this](std::vector<double> &force,
                                                              std::vector<double> &oldForce,
                                                              std::vector<double> &x,
                                                              std::vector<double> &v,
                                                              std::vector<double> &m,
                                                              std::vector<int> &type,
                                                              unsigned long count,
                                                              std::vector<unsigned long> &cell0Items,
                                                              std::vector<unsigned long> &cell1Items,
                                                              std::vector<double> &eps,
                                                              std::vector<double> &sig){
            for(unsigned long indexI : cell0Items){
                for(unsigned long indexJ : cell1Items) {
                    this->fpairFun(force, x, eps, sig, m, type, indexI, indexJ);
                }
            }
        });
    }

    fpair_fun_alt_t FGravityCells::getFastForceAltFunction() {
        return fpairFunAlt;
    }

    fpair_fun_ret_t FGravityCells::getFastForceRetFunction() {
        return fpairFunRet;
    }
} // force