//
// Created by alex on 11.01.23.
//

#include "FLennardJonesCellsOMP.h"
#include "defaults.h"

#include <iostream>

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
        particleContainer.runOnData([&](std::vector<Particle>& particles,
                                        std::vector<Membrane>& membranes,
                                        ParticleContainer::VectorCoordWrapper& cells,
                                        unsigned long count,
                                        std::vector<unsigned long>& activeParticles,
                                        std::unordered_map<unsigned long, unsigned long> &id_to_index){
            using cell_ptr = std::vector<unsigned long> *;
            std::vector<std::vector<cell_ptr>> tasks;
            size_t buffer_size = 0;
            pair_fun_t* pairFun = &this->pairFun;
            cell_ptr c_ptr;
            static const double rt3_2 = std::pow(2, 1 / 3);
            const std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>>& taskGroups = particleContainer.generateDistinctCellNeighbours();
            unsigned long indexI;
            unsigned long indexJ;
            unsigned long indexII;
            unsigned long indexJJ;
            unsigned long indexX;
            unsigned long indexY;
            unsigned long indexZ;
            unsigned long indexC0;
            unsigned long indexC1;

            tasks.emplace_back();
            for (auto &cell: cells) {
                if (buffer_size > max_thread_tasks) {
                    tasks.emplace_back();
                    buffer_size = 0;
                }
                tasks.back().emplace_back(&cell);
                buffer_size += cell.size();
            }

            #pragma omp parallel default(none) shared(particles, tasks, pairFun) private(c_ptr, indexI, indexJ,indexII,indexJJ,indexX,indexY)
            {
                #pragma omp for
                for (indexII = 0; indexII < tasks.size(); indexII++) {
                    for (indexJJ = 0; indexJJ < tasks[indexII].size(); indexJJ++) {
                        c_ptr = tasks[indexII][indexJJ];
                        for (indexX = 0; indexX < c_ptr->size(); indexX++) {
                            for (indexY = indexX + 1; indexY < c_ptr->size(); indexY++) {
                                indexI = (*c_ptr)[indexX];
                                indexJ = (*c_ptr)[indexY];
                                (*pairFun)(particles[indexI], particles[indexJ]);
                            }
                        }
                    }
                }

                #pragma omp barrier
            }

            for(indexX = 0; indexX < 26; indexX++) { //task group index: indexX
                #pragma omp parallel default(none) shared(particles, cells, pairFun, taskGroups) private(indexI, indexJ,indexII,indexJJ,indexX,indexY,indexC0,indexC1,indexZ)
                {
                    #pragma omp for schedule(static,1)
                    for (indexY = 0; indexY < taskGroups[indexX].size(); indexY++) { //task index: indexY - each thread gets one task
                        for (indexZ = 0; indexZ < taskGroups[indexX][indexY].size(); indexZ++) { //pair index: indexII
                            indexC0 = taskGroups[indexX][indexY][indexZ].first;
                            indexC1 = taskGroups[indexX][indexY][indexZ].second;
                            for (indexII = 0; indexII < cells[indexC0].size(); indexII++) {
                                for (indexJJ = 0; indexJJ < cells[indexC1].size(); indexJJ++) {
                                    indexI = cells[indexC0][indexII];
                                    indexJ = cells[indexC1][indexJJ];
                                    (*pairFun)(particles[indexI], particles[indexJ]);
                                }
                            }
                        }
                    }
                    #pragma omp barrier
                }
            }
        });
    }
} // force