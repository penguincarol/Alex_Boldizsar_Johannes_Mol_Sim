//
// Created by alex on 11.01.23.
//

#include "FGravityCells.h"

namespace sim::physics::force {
    void FGravityCells::setPairFun() {
        pairFun = forceDelegate.getForceFunction();
        fpairFun = forceDelegate.getFastForceFunction();
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
        particleContainer.forAllCells([this](std::vector<Particle> &particles,
                                             unsigned long count,
                                             std::vector<unsigned long> &cellItems){
            for(unsigned long indexX = 0; indexX < cellItems.size(); indexX++){
                for(unsigned long indexY = indexX + 1; indexY < cellItems.size(); indexY++) {
                    unsigned long indexI = cellItems[indexX];
                    unsigned long indexJ = cellItems[indexY];
                    this->pairFun(particles[indexI], particles[indexJ]);
                }
            }
        });

        particleContainer.forAllDistinctCellNeighbours([this](std::vector<Particle> &particles,
                                                              unsigned long count,
                                                              std::vector<unsigned long> &cell0Items,
                                                              std::vector<unsigned long> &cell1Items){
            for(unsigned long indexI : cell0Items){
                for(unsigned long indexJ : cell1Items) {
                    this->pairFun(particles[indexI], particles[indexJ]);
                }
            }
        });
    }
} // force