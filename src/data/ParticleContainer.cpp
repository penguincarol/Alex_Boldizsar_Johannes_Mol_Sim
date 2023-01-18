#include "io/output/Logging.h"
#include "ParticleContainer.h"
#include "Particle.h"

#include <vector>
#include <iostream>
#include <numeric>
#include <omp.h>

#pragma region Contructors

ParticleContainer::ParticleContainer() {
    count = 0;
}

ParticleContainer::ParticleContainer(const std::vector<Particle> &buffer) {
    count = buffer.size();
    particles.resize(count);

    //define which particles are still part of the simulation
    activeParticles.resize(count);

    //load particles
    for (unsigned long index{0}; index < count; index++) {
        particles[index] = buffer[index];
        id_to_index[buffer[index].getID()] = index;
        index_to_id[index] = buffer[index].getID();
        activeParticles[index] = buffer[index].getID();
    }
}

ParticleContainer::ParticleContainer(const std::vector<Particle> &buffer, std::array<double, 3> domSize,
                                     double r_cutoff, const std::vector<Membrane>& membranesIn, bool fEOMP) :
        ParticleContainer::ParticleContainer(buffer) {
    eOMP = fEOMP;
    domainSize = domSize;
    x_2_max = domainSize[2];
    x_1_max = domainSize[1];
    x_0_max = domainSize[0];
    x_2_min = 0;
    x_1_min = 0;
    x_0_min = 0;
    membranes = membranesIn;
    //very explicit casts to absolutely make sure, that the rounding is done correctly
    //this implementation uses shorter grids on the side if the numbers are nasty btw
    std::array<double, 3> helperGridDimensions{std::ceil(domSize[0] / r_cutoff), std::ceil(domSize[1] / r_cutoff),
                                               std::ceil(domSize[2] / r_cutoff)};
    gridDimensions = {(unsigned int) helperGridDimensions[0], (unsigned int) helperGridDimensions[1],
                      (unsigned int) helperGridDimensions[2]};

    cells = VectorCoordWrapper(gridDimensions[0]+2, gridDimensions[1]+2, gridDimensions[2]+2);
    this->r_cutoff = (double) r_cutoff;

    if(eOMP) {
        //create padding
        unsigned long newSize = (cells.size() - 1) * padding_count + count;
        particles.resize(newSize);
        //particles are now at the beginning of all vectors
        //by calling update cells will be sorted and moved into cell
    }

    updateCells();

    //halo value
    root6_of_2 = std::pow(2, 1/6);

    if(eOMP){
        initTaskModel();
        initAlternativeTaskModel();}
}

ParticleContainer::ParticleContainer(const std::vector<Particle> &buffer, std::array<double, 2> domainSize,
                                     double r_cutoff, const std::vector<Membrane>& membranesIn, bool fEOMP) :
        ParticleContainer::ParticleContainer(buffer, {domainSize[0], domainSize[1], r_cutoff}, r_cutoff, membranesIn, fEOMP) {};
#pragma endregion

#pragma region Utils

unsigned long ParticleContainer::size() {
    return count;
}

unsigned long ParticleContainer::activeSize() {
    return activeParticles.size();
}

void ParticleContainer::clear() {
    count = 0;
    particles.clear();
    activeParticles.clear();
    cells.clear();
    id_to_index.clear();
    index_to_id.clear();
}

Particle ParticleContainer::getParticle(unsigned long i) {
    return particles[i];
}

std::vector<Membrane>& ParticleContainer::getMembranes() {
    return membranes;
}

std::array<unsigned int, 3> ParticleContainer::getGridDimensions() {
    return gridDimensions;
}

void ParticleContainer::updateCells() {
    // have local copy of current particle indices -> is id_to_index map
    // padding already exists
    // create desired cell indices -> will update cells
    io::output::loggers::general->trace("updateCells called");
    for (auto &cell: cells) cell.clear();
    for (unsigned long id: activeParticles) {
        unsigned long i = id_to_index[id];
        const auto x = particles[i].getX();
        //i am intentionally rounding down with casts from double to unsigned int
        std::array<unsigned int, 3> cellCoordinate = {0,0,0};
        if(x[0] > 0) cellCoordinate[0] = (unsigned int) (x[0] / r_cutoff);
        if(x[1] > 0) cellCoordinate[1] = (unsigned int) (x[1] / r_cutoff);
        if(x[2] > 0) cellCoordinate[2] = (unsigned int) (x[2] / r_cutoff);
        this->cells[cellIndexFromCellCoordinates(cellCoordinate)].emplace_back(id);
    } // now cells contain ID of particle -> need sort particles and replace ID in cell with index

    if(!eOMP) return;

    const unsigned long cellCount = cells.size();
    unsigned long vecIndex = 0;
    for (unsigned long indexCells {0}; indexCells < cellCount; indexCells++) {
        const unsigned long cellItems = cells[indexCells].size();
        for (unsigned long indexC {0}; indexC < cellItems; indexC++) {
            const unsigned long thisID = cells[indexCells][indexC];
            if (index_to_id.contains(vecIndex)) { // need to swap
                const unsigned long otherID = index_to_id[vecIndex];
                swap(thisID, otherID);
            }
            else { //current position is free
                move(id_to_index[thisID], vecIndex, thisID);
            }
            //thisID done -> replace ID in cells with index
            cells[indexCells][indexC] = vecIndex;

            vecIndex++;
        }
        vecIndex += padding_count; // padding
    }
}

void ParticleContainer::swap(unsigned long id0, unsigned long id1) {
    unsigned long index0 = id_to_index[id0];
    unsigned long index1 = id_to_index[id1];
    if(index0==index1) return;

    Particle p = particles[index0];
    particles[index0] = particles[index1];
    particles[index1] = p;

    id_to_index[id0] = index1;
    id_to_index[id1] = index0;
    index_to_id[index0] = id1;
    index_to_id[index1] = id0;
}

void ParticleContainer::move(unsigned long indexSrc, unsigned indexDst, unsigned long id) {
    particles[indexDst] = particles[indexSrc];

    id_to_index[id] = indexDst;
    index_to_id.erase(indexSrc);
    index_to_id[indexDst] = id;
}

void ParticleContainer::deactivateParticles(std::unordered_set<unsigned long> &indices) {
    for(unsigned long ind : indices) {
        id_to_index.erase(index_to_id[ind]);
        index_to_id.erase(ind);
    }
    activeParticles.erase(std::remove_if(activeParticles.begin(), activeParticles.end(), [&](const auto &id) {
        return indices.contains(id_to_index[id]);
    }), activeParticles.end());
}

#pragma endregion

#pragma region Functional



[[maybe_unused]] void ParticleContainer::forAllMembraneSprings(const std::function<void(Particle &p1, Particle &p2, double desiredDistance, double springStrength)> &function){
    for(Membrane& membrane: membranes){
        if(membrane.getMembrNodes().size() <= 0){
            continue;
        }

        std::array<size_t,2> membrDims{membrane.getMembrNodes().size(), membrane.getMembrNodes()[0].size()};

        for(size_t i = 0; i < membrDims[0]; i++){
            for(size_t j = 0; j < membrDims[1] - 1; j++){
                auto id1 = membrane.getMembrNodes()[i][j];
                auto id2 = membrane.getMembrNodes()[i][j+1];

                function(particles[id_to_index[id1]], particles[id_to_index[id2]], membrane.getDesiredDistance(), membrane.getSpringStrength());
            }
        }

        for(size_t i = 0; i < membrDims[0] - 1; i++){
            for(size_t j = 0; j < membrDims[1]; j++){
                auto id1 = membrane.getMembrNodes()[i][j];
                auto id2 = membrane.getMembrNodes()[i+1][j];

                function(particles[id_to_index[id1]], particles[id_to_index[id2]], membrane.getDesiredDistance(), membrane.getSpringStrength());
            }
        }

        for(size_t i = 0; i < membrDims[0]; i++) {
            for (size_t j = 0; j < membrDims[1]; j++) {

                //to top right (thinking about membrane structure)
                if(i+1 < membrDims[0] && j+1 < membrDims[1]){
                    auto id1 = membrane.getMembrNodes()[i][j];
                    auto id2 = membrane.getMembrNodes()[i+1][j+1];

                    function(particles[id_to_index[id1]], particles[id_to_index[id2]], membrane.getDesiredDistance(), membrane.getSpringStrength());
                }

                //to bottom right (thinking about membrane structure)
                if(i+1 < membrDims[0] && j-1 >= 0){
                    auto id1 = membrane.getMembrNodes()[i][j];
                    auto id2 = membrane.getMembrNodes()[i+1][j-1];

                    function(particles[id_to_index[id1]], particles[id_to_index[id2]], membrane.getDesiredDistance(), membrane.getSpringStrength());
                }
            }
        }
    }
}

unsigned int ParticleContainer::cellIndexFromCellCoordinatesFast(unsigned int x0, unsigned int x1, unsigned int x2) {
    return (x0 +
            x1 * gridDimensions[0] +
            x2 * gridDimensions[0] * gridDimensions[1]
    );
}

unsigned int ParticleContainer::cellIndexFromCellCoordinates(std::array<unsigned int, 3> coords) {
//If we decide to change the order of the cells (e.g. by using some fancy 3d space filling curve)
// this method obviously has to be rewritten

//the current version:
//cells[0] corresponds to "cells[0][0][0]"
//cells[1] = "cells[1][0][0]"
//cells[domainSize[0]/r_cutoff] = "cells[0][1][0]"
//cells[(domainSize[0]/r_cutoff)/(domainSize[1]/r_cutoff)] = "cells[0][0][1]"
    unsigned int x0 = std::min(coords[0], gridDimensions[0] - 1);
    unsigned int x1 = std::min(coords[1], gridDimensions[1] - 1);
    unsigned int x2 = std::min(coords[2], gridDimensions[2] - 1);

    return (x0 +
            x1 * gridDimensions[0] +
            x2 * gridDimensions[0] * gridDimensions[1]
    );
}

void ParticleContainer::forAllPairsInSameCell(const std::function<void(Particle &p1, Particle &p2)> &function) {
    for (std::vector<unsigned long> &cellItems: cells) {
        for (unsigned long i = 0; i < cellItems.size(); i++) {
            for (unsigned long j = i + 1; j < cellItems.size(); j++) {
                Particle p1, p2;
                function(particles[i], particles[j]);
            }
        }
    }
}

void ParticleContainer::clearStoreForce() {
    for(unsigned long i {0}; i < count; i++) {
        auto& p = particles[i];
        p.setOldF(p.getF());
        p.setF({0,0,0});
    }
}

void ParticleContainer::initAlternativeTaskModel(){
    //26 TaskGroups (for the 13 cases*2)
    //every taskGroup has the pairs of CellIndices that are independently doable split up into numThreads packages

    //All these variables could be const attributes of class
    const auto numCases = 13;

    using a = std::array<int, 3>;
    constexpr std::array<a, numCases> offsets{a{1,0,0}, a{0,1,0}, a{0,0,1},
                                              a{1,1,0}, a{1,-1,0},
                                              a{1,0,1}, a{1,0,-1}, a{0,1,1}, a{0,1,-1},
                                              a{1,1,1}, a{1,-1,1}, a{1,1,-1}, a{1,-1,-1}};

    auto gD = gridDimensions;
    using b = std::array<unsigned int, 3>;
    const std::array<b, numCases> upperBounds{b{gD[0]-1, gD[1], gD[2]}, b{gD[0], gD[1]-1, gD[2]}, b{gD[0], gD[1], gD[2]-1},
                                              b{gD[0]-1, gD[1]-1, gD[2]}, b{gD[0]-1, gD[1], gD[2]},
                                              b{gD[0]-1, gD[1], gD[2]-1}, b{gD[0]-1, gD[1], gD[2]}, b{gD[0], gD[1]-1, gD[2]-1}, b{gD[0], gD[1]-1, gD[2]},
                                              b{gD[0]-1, gD[1]-1, gD[2]-1}, b{gD[0]-1, gD[1], gD[2]-1}, b{gD[0]-1, gD[1]-1, gD[2]}, b{gD[0]-1, gD[1], gD[2]}};

    constexpr std::array<b, numCases> lowerBounds{b{0,0,0}, b{0,0,0}, b{0,0,0},
                                                  b{0,0,0}, b{0,1,0},
                                                  b{0,0,0}, b{0,0,1}, b{0,0,0}, b{0,0,1},
                                                  b{0,0,0}, b{0,1,0}, b{0,0,1}, b{0,1,1}};

    for(auto c = 0; c < numCases; c++){ //pun intended
        //const std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>>& generateDistinctCellNeighbours()
        const unsigned long maxThreads{static_cast<unsigned long>(omp_get_max_threads())};
        std::vector<std::pair<unsigned long, unsigned long>> independentTasksBlock{maxThreads};

        for(unsigned int x0 = lowerBounds[c][0]; x0 < upperBounds[c][0]; x0+=2){
            for(unsigned int x1 = lowerBounds[c][1]; x1 < upperBounds[c][1]; x1+=2){
                for(unsigned int x2 = lowerBounds[c][2]; x2 < upperBounds[c][2]; x2+=2){
                    alternativeTaskModelCache.emplace_back(cellIndexFromCellCoordinatesFast(x0, x1, x2),
                                                       cellIndexFromCellCoordinatesFast(x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2]));
                    SPDLOG_TRACE("Added CellInteraction (({} {} {}), ({} {} {})) to taskBlock {}", x0, x1, x2, x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2], 2*c+0);
                }
            }
        }

        //yes you could cut this down to 2 lines with another helper array but this more verbose version seems much easier to understand
        std::vector<std::pair<unsigned long, unsigned long>> independentTasksBlock2{maxThreads};

        for(unsigned int x0 = lowerBounds[c][0] + 1; x0 < upperBounds[c][0]; x0+=2){
            for(unsigned int x1 = lowerBounds[c][1] + 1; x1 < upperBounds[c][1]; x1+=2){
                for(unsigned int x2 = lowerBounds[c][2] + 1; x2 < upperBounds[c][2]; x2+=2){
                    alternativeTaskModelCache.emplace_back(cellIndexFromCellCoordinatesFast(x0, x1, x2),
                                                       cellIndexFromCellCoordinatesFast(x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2]));
                    SPDLOG_TRACE("Added CellInteraction (({} {} {}), ({} {} {})) to taskBlock {}", x0, x1, x2, x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2], 2*c+1);
                }
            }
        }


    }
}

void ParticleContainer::initTaskModel() {
    //26 TaskGroups (for the 13 cases*2)
    //every taskGroup has the pairs of CellIndices that are independently doable split up into numThreads packages

    //All these variables could be const attributes of class
    const auto numCases = 13;

    using a = std::array<int, 3>;
    constexpr std::array<a, numCases> offsets{a{1,0,0}, a{0,1,0}, a{0,0,1},
                                              a{1,1,0}, a{1,-1,0},
                                              a{1,0,1}, a{1,0,-1}, a{0,1,1}, a{0,1,-1},
                                              a{1,1,1}, a{1,-1,1}, a{1,1,-1}, a{1,-1,-1}};

    auto gD = gridDimensions;
    using b = std::array<unsigned int, 3>;
    const std::array<b, numCases> upperBounds{b{gD[0]-1, gD[1], gD[2]}, b{gD[0], gD[1]-1, gD[2]}, b{gD[0], gD[1], gD[2]-1},
                                              b{gD[0]-1, gD[1]-1, gD[2]}, b{gD[0]-1, gD[1], gD[2]},
                                              b{gD[0]-1, gD[1], gD[2]-1}, b{gD[0]-1, gD[1], gD[2]}, b{gD[0], gD[1]-1, gD[2]-1}, b{gD[0], gD[1]-1, gD[2]},
                                              b{gD[0]-1, gD[1]-1, gD[2]-1}, b{gD[0]-1, gD[1], gD[2]-1}, b{gD[0]-1, gD[1]-1, gD[2]}, b{gD[0]-1, gD[1], gD[2]}};

    constexpr std::array<b, numCases> lowerBounds{b{0,0,0}, b{0,0,0}, b{0,0,0},
                                                  b{0,0,0}, b{0,1,0},
                                                  b{0,0,0}, b{0,0,1}, b{0,0,0}, b{0,0,1},
                                                  b{0,0,0}, b{0,1,0}, b{0,0,1}, b{0,1,1}};

    for(auto c = 0; c < numCases; c++){ //pun intended
        //const std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>>& generateDistinctCellNeighbours()
        const unsigned long maxThreads{static_cast<unsigned long>(omp_get_max_threads())};
        std::vector<std::vector<std::pair<unsigned long, unsigned long>>> independentTasksBlock{maxThreads, std::vector<std::pair<unsigned long, unsigned long>>{}};

        size_t roundRobinCounter{0};
        for(unsigned int x0 = lowerBounds[c][0]; x0 < upperBounds[c][0]; x0+=2){
            for(unsigned int x1 = lowerBounds[c][1]; x1 < upperBounds[c][1]; x1+=2){
                for(unsigned int x2 = lowerBounds[c][2]; x2 < upperBounds[c][2]; x2+=2){
                    independentTasksBlock[roundRobinCounter].emplace_back(cellIndexFromCellCoordinatesFast(x0, x1, x2),
                                                                          cellIndexFromCellCoordinatesFast(x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2]));
                    roundRobinCounter = (roundRobinCounter+1)%maxThreads;
                    SPDLOG_TRACE("Added CellInteraction (({} {} {}), ({} {} {})) to taskBlock {}", x0, x1, x2, x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2], 2*c+0);
                }
            }
        }

        taskModelCache.emplace_back(independentTasksBlock);

        //yes you could cut this down to 2 lines with another helper array but this more verbose version seems much easier to understand
        std::vector<std::vector<std::pair<unsigned long, unsigned long>>> independentTasksBlock2{maxThreads, std::vector<std::pair<unsigned long, unsigned long>>{}};
        roundRobinCounter = 0;
        for(unsigned int x0 = lowerBounds[c][0] + 1; x0 < upperBounds[c][0]; x0+=2){
            for(unsigned int x1 = lowerBounds[c][1] + 1; x1 < upperBounds[c][1]; x1+=2){
                for(unsigned int x2 = lowerBounds[c][2] + 1; x2 < upperBounds[c][2]; x2+=2){
                    independentTasksBlock2[roundRobinCounter].emplace_back(cellIndexFromCellCoordinatesFast(x0, x1, x2),
                                                                           cellIndexFromCellCoordinatesFast(x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2]));
                    roundRobinCounter = (roundRobinCounter+1)%maxThreads;
                    SPDLOG_TRACE("Added CellInteraction (({} {} {}), ({} {} {})) to taskBlock {}", x0, x1, x2, x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2], 2*c+1);
                }
            }
        }
        taskModelCache.emplace_back(independentTasksBlock2);

    }
}

const std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>>& ParticleContainer::generateDistinctCellNeighbours() {
    return taskModelCache;
}

const std::vector<std::pair<unsigned long, unsigned long>>& ParticleContainer::generateDistinctAlternativeCellNeighbours(){
    return alternativeTaskModelCache;
}

#pragma endregion

