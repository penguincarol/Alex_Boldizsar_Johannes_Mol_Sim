//
// Created by alex on 11.01.23.
//

#include "FLennardJonesCellsOMP.h"
#include "defaults.h"

#include <iostream>

template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &orig)
{
    std::vector<T> ret;
    for(const auto &v: orig)
        ret.insert(ret.end(), v.begin(), v.end());
    return ret;
}

namespace sim::physics::force {
    /**
     * @brief Calculates the using the Linked-Cell algorithm. The linked cell algorithm has a much better runtime than the All-Pairs algorithm (see plot)
     * \image html plot.png width=800px
     * \image latex plot.eps "Runtime comparison of All-Pairs algorithm with Linked-Cell algorithm" width=5cm
     */
    void FLennardJonesCellsOMP::operator()() {
        /*
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
        });*/
        particleContainer.runOnDataCell([&](std::vector<double> &force,
                                            std::vector<double> &oldForce,
                                            std::vector<double> &x,
                                            std::vector<double> &v,
                                            std::vector<double> &m,
                                            std::vector<int> &type,
                                            unsigned long count,
                                            ParticleContainer::VectorCoordWrapper& cells,
                                            std::vector<double> &eps,
                                            std::vector<double> &sig){
            //size_t interactions{0};
            #pragma omp parallel for //reduction(+:interactions)
            for(size_t cellIndex=0; cellIndex < cells.size(); cellIndex++){
                auto& cell = cells[cellIndex];
                for(size_t i = 0; i < cell.size(); i++){
                    for(size_t j = i+1; j < cell.size(); j++){
                        //interactions++;
                        this->fpairFun(force, x, eps, sig, m, type, i, j);
                    }
                }
            }
            #pragma omp barrier
            //std::cout<<"Inter Cells with themselves: " << interactions << std::endl;
        });

        //size_t interactionsDistinctCells{0};

        std::vector<std::vector<std::pair<unsigned long, unsigned long>>> taskBlocks = particleContainer.generateDistinctAlternativeCellNeighbours();
        std::vector<std::pair<unsigned long, unsigned long>> tasks = flatten(taskBlocks);

            particleContainer.runOnDataCell([&tasks, this](std::vector<double> &force,
                                                std::vector<double> &oldForce,
                                                std::vector<double> &x,
                                                std::vector<double> &v,
                                                std::vector<double> &m,
                                                std::vector<int> &type,
                                                unsigned long count,
                                                ParticleContainer::VectorCoordWrapper& cells,
                                                std::vector<double> &eps,
                                                std::vector<double> &sig){
            #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

            #pragma omp parallel for reduction(vec_double_plus : force)
                for(auto& [cellIndexI, cellIndexJ]: tasks) {
                    size_t cII = cellIndexI;     //openMP problems with the other variant
                    size_t cIJ = cellIndexJ;

                    for (auto pIndexI: cells[cII]) {
                        for (auto pIndexJ: cells[cIJ]) {
                            this->fpairFun(force, x, eps, sig, m, type, pIndexI, pIndexJ);
                            //interactionsDistinctCells++;
                        }
                    }
                }
            });

        /*
        size_t interactionsDistinctCells_st{0};
        particleContainer.forAllDistinctCellNeighbours([this, &interactionsDistinctCells_st](std::vector<double> &force,
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
                    //this->fpairFun(force, x, eps, sig, m, type, indexI, indexJ);
                    interactionsDistinctCells_st++;
                }
            }
        });*/

        //std::cout<<"Inter between cells old version: " << interactionsDistinctCells_st << " new version: "<< interactionsDistinctCells << std::endl;
    }

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
    }

    fpair_fun_t FLennardJonesCellsOMP::getFastForceFunction() {
        return fpairFun;
    }
} // sim::physics::force
