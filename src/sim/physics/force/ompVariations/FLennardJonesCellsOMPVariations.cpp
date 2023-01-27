#include "sim/physics/force/FLennardJonesCellsOMP.h"
#include "defaults.h"

#include <iostream>

namespace sim::physics::force {

#ifdef ONE_DIMENSIONAL_TASKS
    /**
 * I am not writing yet another task initializer so this method just creates the 1d tasks using the 2d tasks
 */
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>> &orig)
{
    std::vector<T> ret;
    for(const auto &v: orig)
        ret.insert(ret.end(), v.begin(), v.end());
    return ret;
}


    /**
     * This is the implementaiton of the operator utilizing the one-dimensional-task approach as displayed in the presentation
     */
    void FLennardJonesCellsOMP::operator()() {

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
            auto fPairFun = this->fpairFun;
            #pragma omp parallel for default(none) shared(fPairFun, force, oldForce, x, v, m, type, count, cells, eps, sig)//reduction(+:interactions)
            for(size_t cellIndex=0; cellIndex < cells.size(); cellIndex++){
                auto& cell = cells[cellIndex];
                for(size_t i = 0; i < cell.size(); i++){
                    for(size_t j = i+1; j < cell.size(); j++){
                        //interactions++;
                        fPairFun(force, x, eps, sig, m, type, i, j);
                    }
                }
            }
            //std::cout<<"Inter Cells with themselves: " << interactions << std::endl;
        });

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
            auto fPairFun = this->fpairFun;
            #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

            #pragma omp parallel for default(none) shared(fPairFun, oldForce, x, v, m, type, count, cells, eps, sig, tasks) reduction(vec_double_plus : force)
            for(auto& [cellIndexI, cellIndexJ]: tasks) {
                size_t cII = cellIndexI;     //openMP problems with the other variant
                size_t cIJ = cellIndexJ;

                for (auto pIndexI: cells[cII]) {
                    for (auto pIndexJ: cells[cIJ]) {
                        fPairFun(force, x, eps, sig, m, type, pIndexI, pIndexJ);
                        //interactionsDistinctCells++;
                    }
                }
            }
        });
    }

#endif

#ifdef THREE_DIMENSIONAL_TASKS

    /**
     * This is the implementaiton of the operator utilizing the three-dimensional-task approach as displayed in the presentation
     */
    void FLennardJonesCellsOMP::operator()() {
        particleContainer.runOnDataCell([&](std::vector<double> &force,
                                            std::vector<double> &oldForce,
                                            std::vector<double> &x,
                                            std::vector<double> &v,
                                            std::vector<double> &m,
                                            std::vector<int> &type,
                                            unsigned long count,
                                            ParticleContainer::VectorCoordWrapper &cells,
                                            std::vector<double> &eps,
                                            std::vector<double> &sig) {

        auto fPairFun = this->fpairFun;
            //size_t interactions{0};
#pragma omp parallel for default(none) shared(fPairFun, force, oldForce, x, v, m, type, count, cells, eps, sig)//reduction(+:interactions)
            for (size_t cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
                auto &cell = cells[cellIndex];
                for (size_t i = 0; i < cell.size(); i++) {
                    for (size_t j = i + 1; j < cell.size(); j++) {
                        //interactions++;
                        fPairFun(force, x, eps, sig, m, type, i, j);
                    }
                }
            }
            //std::cout<<"Inter Cells with themselves: " << interactions << std::endl;
        });

        std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>> taskBlocks = particleContainer.generateDistinctCellNeighbours();
        auto fPairFun = this->fpairFun;
        for (auto &tasks: taskBlocks) {
            //the previous taskBlock needs to be finished before the next taskBlock can start>
            if (tasks.size() != omp_get_max_threads()) {
                io::output::loggers::simulation->debug(
                        "Task creation for force calculation didn't result in appropriate amount of tasks");
            }
#           pragma omp parallel for default(none) shared(tasks, fPairFun) schedule(static, 1)
            for (auto &task: tasks) {
                for (auto &[cellIndexI, cellIndexJ]: task) {
                    size_t cII = cellIndexI;     //openMP problems with the other variant
                    size_t cIJ = cellIndexJ;

                    particleContainer.runOnDataCell([&](std::vector<double> &force,
                                                        std::vector<double> &oldForce,
                                                        std::vector<double> &x,
                                                        std::vector<double> &v,
                                                        std::vector<double> &m,
                                                        std::vector<int> &type,
                                                        unsigned long count,
                                                        ParticleContainer::VectorCoordWrapper &cells,
                                                        std::vector<double> &eps,
                                                        std::vector<double> &sig) {
                        for (auto pIndexI: cells[cII]) {
                            for (auto pIndexJ: cells[cIJ]) {
                                fPairFun(force, x, eps, sig, m, type, pIndexI, pIndexJ);
                            }
                        }
                    });
                }
            }
        }
        //#pragma omp barrier
    }
#endif
}