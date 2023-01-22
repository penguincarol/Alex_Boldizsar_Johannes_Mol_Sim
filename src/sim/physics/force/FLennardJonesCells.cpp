#include "FLennardJonesCells.h"
#include "FLennardJones.h"

namespace sim::physics::force {
    /**
     * @brief Calculates the using the Linked-Cell algorithm. The linked cell algorithm has a much better runtime than the All-Pairs algorithm (see plot)
     * \image html plot.png width=800px
     * \image latex plot.eps "Runtime comparison of All-Pairs algorithm with Linked-Cell algorithm" width=5cm
     */
    void FLennardJonesCells::operator()() {
        particleContainer.forAllCells([this](vec4d_t &force,
                                             vec4d_t &oldForce,
                                             vec4d_t &x,
                                             vec4d_t &v,
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

        particleContainer.forAllDistinctCellNeighbours([this](vec4d_t &force,
                                                              vec4d_t &oldForce,
                                                              vec4d_t &x,
                                                              vec4d_t &v,
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
        fpairFunAlt = forceDelegate.getFastForceAltFunction();
        fpairFunRet = forceDelegate.getFastForceRetFunction();
    }

    fpair_fun_ret_t FLennardJonesCells::getFastForceRetFunction() {
        return fpairFunRet;
    }

    fpair_fun_t FLennardJonesCells::getFastForceFunction() {
        return fpairFun;
    }

    fpair_fun_alt_t FLennardJonesCells::getFastForceAltFunction() {
        return fpairFunAlt;
    }
} // sim::physics::force
