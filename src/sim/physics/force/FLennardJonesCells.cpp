#include "FLennardJonesCells.h"
#include "FLennardJones.h"

namespace sim::physics::force {
    /**
     * @brief Calculates the using the Linked-Cell algorithm. The linked cell algorithm has a much better runtime than the All-Pairs algorithm (see plot)
     * \image html plot.png width=800px
     * \image latex plot.eps "Runtime comparison of All-Pairs algorithm with Linked-Cell algorithm" width=5cm
     */
    void FLennardJonesCells::operator()() {
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

    /**
     * @brief Returns the force function used
     * 
     * @return pair_fun_t& 
     */
    pair_fun_t &FLennardJonesCells::getForceFunction() {
        return pairFun;
    }

    /**
     * @brief The name says it all
     * 
     * @param pc 
     */
    void FLennardJonesCells::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
        forceDelegate.setParticleContainer(pc);
        setPairFun();
    }

    void FLennardJonesCells::setPairFun() {
        pairFun = forceDelegate.getForceFunction();
        fpairFun = forceDelegate.getFastForceFunction();
    }

    fpair_fun_t FLennardJonesCells::getFastForceFunction() {
        return fpairFun;
    }
} // sim::physics::force
