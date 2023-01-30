#include "sim/physics/force/FLennardJonesCellsOMP.h"
#include "defaults.h"

#include <iostream>
#include <omp.h>

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
        //io::output::loggers::general->error("Jup, one-dim used!");
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
                        fPairFun(force, x, eps, sig, m, type, cell[i], cell[j]);
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
            static const double rt3_2 = std::pow(2, 1.0 / 3.0);
            auto fPairFun = this->fpairFun;
            #pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

            #pragma omp parallel for default(none) shared(rt3_2, fPairFun, oldForce, x, v, m, type, count, cells, eps, sig, tasks) reduction(vec_double_plus : force)
            for(auto& [cellIndexI, cellIndexJ]: tasks) {
                size_t cII = cellIndexI;     //openMP problems with the other variant
                size_t cIJ = cellIndexJ;

                for (auto pIndexI: cells[cII]) {
                    for (auto pIndexJ: cells[cIJ]) {

                        auto dsqr = ((x[3*pIndexI + 0]-x[3*pIndexJ+0])*(x[3*pIndexI + 0]-x[3*pIndexJ+0]))+
                                ((x[3*pIndexI + 1]-x[3*pIndexJ+1])*(x[3*pIndexI + 1]-x[3*pIndexJ+1]))+
                                ((x[3*pIndexI + 2]-x[3*pIndexJ+2])*(x[3*pIndexI + 2]-x[3*pIndexJ+2]));
                        auto sigma2 = (sig[pIndexI] + sig[pIndexJ])*(sig[pIndexI] + sig[pIndexJ]) / 4;  //sigma^2
                        if (type[pIndexI] & 0x80000000 || type[pIndexJ] & 0x80000000) {
                            if (dsqr >= rt3_2 * sigma2) continue;
                        }
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
        //io::output::loggers::general->error("Jup, three-dim used!");
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
                        fPairFun(force, x, eps, sig, m, type, cell[i], cell[j]);
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
                        static const double rt3_2 = std::pow(2, 1.0 / 3.0);
                        for (auto pIndexI: cells[cII]) {
                            for (auto pIndexJ: cells[cIJ]) {

                                auto dsqr = ((x[3*pIndexI + 0]-x[3*pIndexJ+0])*(x[3*pIndexI + 0]-x[3*pIndexJ+0]))+
                                            ((x[3*pIndexI + 1]-x[3*pIndexJ+1])*(x[3*pIndexI + 1]-x[3*pIndexJ+1]))+
                                            ((x[3*pIndexI + 2]-x[3*pIndexJ+2])*(x[3*pIndexI + 2]-x[3*pIndexJ+2]));
                                auto sigma2 = (sig[pIndexI] + sig[pIndexJ])*(sig[pIndexI] + sig[pIndexJ]) / 4;  //sigma^2
                                if (type[pIndexI] & 0x80000000 || type[pIndexJ] & 0x80000000) {
                                    if (dsqr >= rt3_2 * sigma2) continue;
                                }

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

#ifdef TASK_ORIENTED_2D
    void FLennardJonesCellsOMP::operator()() {
        //io::output::loggers::general->error("Jup, task oriented 2D is used!");
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
            auto fpairFun = this->fpairFun;
            #pragma omp parallel for default(none) shared(cells, x, eps, sig, m, type, force, fpairFun)
            for (size_t cellIndex = 0; cellIndex < cells.size(); cellIndex++) {
                auto &cell = cells[cellIndex];
                for (size_t i = 0; i < cell.size(); i++) {
                    for (size_t j = i + 1; j < cell.size(); j++) {
                        fpairFun(force, x, eps, sig, m, type, cell[i], cell[j]);
                    }
                }
            }
        });

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

            const std::vector<std::vector<std::pair<unsigned long, unsigned long>>> &taskOrientedGroups = particleContainer.generateDistinctTaskOrientedCellNeighbours();
            for(auto& tasks : taskOrientedGroups){
                auto fpairFun = this->fpairFun;
            #pragma omp parallel for default(none) shared(fpairFun, force, oldForce, x, v, m, type, count, cells, eps, sig, tasks)
                for(auto &[cellIndexI, cellIndexJ]: tasks){
                    auto& cellI = cells[cellIndexI];
                    auto& cellJ = cells[cellIndexJ];

                    for(auto& pIndexI: cellI){
                        for(auto& pIndexJ: cellJ){
                            fpairFun(force, x, eps, sig, m, type, pIndexI, pIndexJ);
                        }
                    }

                }
            }
        });

    }
#endif
}

