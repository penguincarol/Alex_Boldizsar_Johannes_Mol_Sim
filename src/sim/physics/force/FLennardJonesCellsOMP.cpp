//
// Created by alex on 11.01.23.
//

#include "FLennardJonesCellsOMP.h"
#include "defaults.h"

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
    }

    fpair_fun_t FLennardJonesCellsOMP::getFastForceFunction() {
        return fpairFun;
    }

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
            //generate tasks: for each individual cell
            using cell_ptr = std::vector<unsigned long>*;
            std::vector<std::vector<cell_ptr>> tasks;
            tasks.emplace_back();
            size_t buffer_size = 0;
            for (auto& cell : cells) {
                if(buffer_size > max_thread_tasks) {
                    tasks.emplace_back();
                    buffer_size = 0;
                }
                tasks.back().emplace_back(&cell);
                buffer_size += cell.size();
            }
            auto fpairFun = this->fpairFun;
            #pragma omp parallel for default(none) shared(tasks, force, x, eps, sig, m, type, fpairFun)
            for (auto& task : tasks) { //task is vector of cell ptr
                for (cell_ptr c_ptr : task) {
                    for(unsigned long indexX = 0; indexX < c_ptr->size(); indexX++){
                        for(unsigned long indexY = indexX + 1; indexY < c_ptr->size(); indexY++) {
                            unsigned long indexI = (*c_ptr)[indexX];
                            unsigned long indexJ = (*c_ptr)[indexY];
                            fpairFun(force, x, eps, sig, m, type, indexI, indexJ);
                        }
                    }
                }
            }

            //generate tasks: for all distinct cell neighbours
            std::vector<std::vector<std::pair<unsigned long, unsigned long>>> alternativeTaskGroups = particleContainer.generateDistinctAlternativeCellNeighbours();
            for(auto& task : alternativeTaskGroups){
                #pragma omp parallel for reduction (+:force) readonly(alternativeTaskGroups, cells, force, x, eps, sig, m, type, fpairFun)

                for(auto [indexC0, indexC1] : task) {
                    auto& cell0 = cells[indexC0];
                    auto& cell1 = cells[indexC1];

                    for(unsigned long indexI : cell0){
                        for(unsigned long indexJ : cell1) {
                            fpairFun(force, x, eps, sig, m, type, indexI, indexJ);
                        }
                    }
                }

            }

            /*
            std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>> task_groups =
                    particleContainer.generateDistinctCellNeighbours();
            for(auto& task_group : task_groups) {
                #pragma omp parallel for default(none) shared(task_group, cells, force, x, eps, sig, m, type, fpairFun)
                for(auto& task : task_group) {
                    for(auto [indexC0, indexC1] : task) {
                        auto& cell0 = cells[indexC0];
                        auto& cell1 = cells[indexC1];

                        for(unsigned long indexI : cell0){
                            for(unsigned long indexJ : cell1) {
                                fpairFun(force, x, eps, sig, m, type, indexI, indexJ);
                            }
                        }
                    }
                }
            }*/
        });
    }
} // force