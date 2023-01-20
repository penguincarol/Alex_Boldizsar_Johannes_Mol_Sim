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
    force.resize(count * 3);
    oldForce.resize(count * 3);
    x.resize(count * 3);
    v.resize(count * 3);
    m.resize(count);
    type.resize(count);
    eps.resize(count);
    sig.resize(count);

    //define which particles are still part of the simulation
    activeParticles.resize(count);

    //load particles
    for (unsigned long index{0}; index < count; index++) {
        auto &f = buffer[index].getF();
        force[index * 3 + 0] = f[0];
        force[index * 3 + 1] = f[1];
        force[index * 3 + 2] = f[2];

        auto &of = buffer[index].getOldF();
        oldForce[index * 3 + 0] = of[0];
        oldForce[index * 3 + 1] = of[1];
        oldForce[index * 3 + 2] = of[2];

        auto &xx = buffer[index].getX();
        x[index * 3 + 0] = xx[0];
        x[index * 3 + 1] = xx[1];
        x[index * 3 + 2] = xx[2];

        auto &vv = buffer[index].getV();
        v[index * 3 + 0] = vv[0];
        v[index * 3 + 1] = vv[1];
        v[index * 3 + 2] = vv[2];

        m[index] = buffer[index].getM();
        type[index] = buffer[index].getType();

        sig[index] = buffer[index].getSigma();
        eps[index] = buffer[index].getEpsilon();

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
        force.resize(newSize * 3);
        oldForce.resize(newSize * 3);
        x.resize(newSize * 3);
        v.resize(newSize * 3);
        m.resize(newSize);
        type.resize(newSize);
        eps.resize(newSize);
        sig.resize(newSize);
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
    force.clear();
    oldForce.clear();
    x.clear();
    v.clear();
    m.clear();
    type.clear();
    eps.clear();
    sig.clear();
    activeParticles.clear();
    cells.clear();
    id_to_index.clear();
    index_to_id.clear();
}

Particle ParticleContainer::getParticle(unsigned long i) {
    Particle p;
    loadParticle(p, i);
    return p;
}

std::vector<Membrane>& ParticleContainer::getMembranes() {
    return membranes;
}

std::array<unsigned int, 3> ParticleContainer::getGridDimensions() {
    return gridDimensions;
}

void ParticleContainer::loadParticle(Particle &p, unsigned long index, std::vector<double> &force,
                                     std::vector<double> &oldForce, std::vector<double> &x, std::vector<double> &v,
                                     std::vector<double> &m, std::vector<int> &type, std::vector<double>& e,
                                     std::vector<double> &s) {
    Eigen::Vector3d f{force[index * 3 + 0],
                      force[index * 3 + 1],
                      force[index * 3 + 2]};
    p.setF(f);
    Eigen::Vector3d of{oldForce[index * 3 + 0],
                       oldForce[index * 3 + 1],
                       oldForce[index * 3 + 2]};
    p.setOldF(of);
    Eigen::Vector3d xx{x[index * 3 + 0],
                       x[index * 3 + 1],
                       x[index * 3 + 2]};
    p.setX(xx);
    Eigen::Vector3d vv{v[index * 3 + 0],
                       v[index * 3 + 1],
                       v[index * 3 + 2]};
    p.setV(vv);
    p.setM(m[index]);
    p.setType(type[index]);
    p.setEpsilon(e[index]);
    p.setSigma(s[index]);
}

void ParticleContainer::loadParticle(Particle &p, unsigned long index) {
    loadParticle(p, index, force, oldForce, x, v, m, type, eps, sig);
}

void ParticleContainer::storeParticle(Particle &p, unsigned long index, std::vector<double> &force,
                                      std::vector<double> &oldForce, std::vector<double> &x, std::vector<double> &v,
                                      std::vector<double> &m, std::vector<int> &type, std::vector<double>& e,
                                      std::vector<double> &s) {
    auto &ff = p.getF();
    force[index * 3 + 0] = ff[0];
    force[index * 3 + 1] = ff[1];
    force[index * 3 + 2] = ff[2];

    auto &oof = p.getOldF();
    oldForce[index * 3 + 0] = oof[0];
    oldForce[index * 3 + 1] = oof[1];
    oldForce[index * 3 + 2] = oof[2];

    auto &xxx = p.getX();
    x[index * 3 + 0] = xxx[0];
    x[index * 3 + 1] = xxx[1];
    x[index * 3 + 2] = xxx[2];

    auto &vvv = p.getV();
    v[index * 3 + 0] = vvv[0];
    v[index * 3 + 1] = vvv[1];
    v[index * 3 + 2] = vvv[2];

    m[index] = p.getM();
    type[index] = p.getType();
    e[index] = p.getEpsilon();
    s[index] = p.getSigma();

}

void ParticleContainer::storeParticle(Particle &p, unsigned long index) {
    storeParticle(p, index, force, oldForce, x, v, m, type, eps, sig);
}


void ParticleContainer::updateCells() {
    // have local copy of current particle indices -> is id_to_index map
    // padding already exists
    // create desired cell indices -> will update cells
    io::output::loggers::general->trace("updateCells called");
    for (auto &cell: cells) cell.clear();
    for (unsigned long id: activeParticles) {
        unsigned long i = id_to_index[id];
        //i am intentionally rounding down with casts from double to unsigned int
        std::array<unsigned int, 3> cellCoordinate = {0,0,0};
        if(x[3*i+0] > 0) cellCoordinate[0] = (unsigned int) (x[3 * i] / r_cutoff);
        if(x[3*i+1] > 0) cellCoordinate[1] = (unsigned int) (x[3 * i+1] / r_cutoff);
        if(x[3*i+2] > 0) cellCoordinate[2] = (unsigned int) (x[3 * i+2] / r_cutoff);
        this->cells[cellIndexFromCellCoordinates(cellCoordinate)].emplace_back(id);
    } // now cells contain ID of particle -> need sort particles and replace ID in cell with index

    if(!false) return;

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

    double f0, f1, f2, of0, of1, of2, x0, x1, x2, v0, v1, v2, mm, s, e;
    int t;

    f0 = force[index0 * 3 + 0];
    f1 = force[index0 * 3 + 1];
    f2 = force[index0 * 3 + 2];
    of0 = oldForce[index0 * 3 + 0];
    of1 = oldForce[index0 * 3 + 1];
    of2 = oldForce[index0 * 3 + 2];
    x0 = x[index0 * 3 + 0];
    x1 = x[index0 * 3 + 1];
    x2 = x[index0 * 3 + 2];
    v0 = v[index0 * 3 + 0];
    v1 = v[index0 * 3 + 1];
    v2 = v[index0 * 3 + 2];
    mm = m[index0];
    t = type[index0];
    s = sig[index0];
    e = eps[index0];

    force[index0 * 3 + 0]    = force[index1 * 3 + 0];
    force[index0 * 3 + 1]    = force[index1 * 3 + 1];
    force[index0 * 3 + 2]    = force[index1 * 3 + 2];
    oldForce[index0 * 3 + 0] = oldForce[index1 * 3 + 0];
    oldForce[index0 * 3 + 1] = oldForce[index1 * 3 + 1];
    oldForce[index0 * 3 + 2] = oldForce[index1 * 3 + 2];
    x[index0 * 3 + 0]        = x[index1 * 3 + 0];
    x[index0 * 3 + 1]        = x[index1 * 3 + 1];
    x[index0 * 3 + 2]        = x[index1 * 3 + 2];
    v[index0 * 3 + 0]        = v[index1 * 3 + 0];
    v[index0 * 3 + 1]        = v[index1 * 3 + 1];
    v[index0 * 3 + 2]        = v[index1 * 3 + 2];
    m[index0]                = m[index1];
    type[index0]             = type[index1];
    sig[index0]              = sig[index1];
    eps[index0]              = eps[index1];

    force[index1 * 3 + 0] = f0;
    force[index1 * 3 + 1] = f1;
    force[index1 * 3 + 2] = f2;
    oldForce[index1 * 3 + 0] = of0;
    oldForce[index1 * 3 + 1] = of1;
    oldForce[index1 * 3 + 2] = of2;
    x[index1 * 3 + 0] = x0;
    x[index1 * 3 + 1] = x1;
    x[index1 * 3 + 2] = x2;
    v[index1 * 3 + 0] = v0;
    v[index1 * 3 + 1] = v1;
    v[index1 * 3 + 2] = v2;
    m[index1] = mm;
    type[index1] = t;
    sig[index1] = s;
    eps[index1] = e;

    id_to_index[id0] = index1;
    id_to_index[id1] = index0;
    index_to_id[index0] = id1;
    index_to_id[index1] = id0;
}

void ParticleContainer::move(unsigned long indexSrc, unsigned indexDst, unsigned long id) {
    force[indexDst * 3 + 0]    = force[indexSrc * 3 + 0];
    force[indexDst * 3 + 1]    = force[indexSrc * 3 + 1];
    force[indexDst * 3 + 2]    = force[indexSrc * 3 + 2];
    oldForce[indexDst * 3 + 0] = oldForce[indexSrc * 3 + 0];
    oldForce[indexDst * 3 + 1] = oldForce[indexSrc * 3 + 1];
    oldForce[indexDst * 3 + 2] = oldForce[indexSrc * 3 + 2];
    x[indexDst * 3 + 0]        = x[indexSrc * 3 + 0];
    x[indexDst * 3 + 1]        = x[indexSrc * 3 + 1];
    x[indexDst * 3 + 2]        = x[indexSrc * 3 + 2];
    v[indexDst * 3 + 0]        = v[indexSrc * 3 + 0];
    v[indexDst * 3 + 1]        = v[indexSrc * 3 + 1];
    v[indexDst * 3 + 2]        = v[indexSrc * 3 + 2];
    m[indexDst]                = m[indexSrc];
    type[indexDst]             = type[indexSrc];
    sig[indexDst]              = sig[indexSrc];
    eps[indexDst]              = eps[indexSrc];

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

void ParticleContainer::forAllParticles(const std::function<void(Particle &)> &function) {
    for (unsigned long id: activeParticles) {
        Particle p;
        loadParticle(p, id_to_index[id]);
        function(p);
        storeParticle(p, id_to_index[id]);
    }
}

void ParticleContainer::forAllParticles(void(*function)(Particle &)) {
    for (unsigned long id: activeParticles) {
        Particle p;
        loadParticle(p, id_to_index[id]);
        function(p);
        storeParticle(p, id_to_index[id]);
    }
}

void ParticleContainer::forAllPairs(void (*function)(Particle &p1, Particle &p2)) {
    for (u_int32_t i = 0; i < activeSize(); i++) {
        for (u_int32_t j = i + 1; j < activeSize(); j++) {
            Particle p1;
            loadParticle(p1, id_to_index[activeParticles[i]]);
            Particle p2;
            loadParticle(p2, id_to_index[activeParticles[j]]);
            function(p1, p2);
            storeParticle(p1, id_to_index[activeParticles[i]]);
            storeParticle(p2, id_to_index[activeParticles[j]]);
        }
    }
}

void ParticleContainer::forAllPairs(const std::function<void(Particle &p1, Particle &p2)> &function) {
    for (u_int32_t i = 0; i < activeSize(); i++) {
        for (u_int32_t j = i + 1; j < activeSize(); j++) {
            Particle p1;
            loadParticle(p1, id_to_index[activeParticles[i]]);
            Particle p2;
            loadParticle(p2, id_to_index[activeParticles[j]]);
            function(p1, p2);
            storeParticle(p1, id_to_index[activeParticles[i]]);
            storeParticle(p2, id_to_index[activeParticles[j]]);
        }
    }
}


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

                Particle p1;
                loadParticle(p1, id_to_index[id1]);
                Particle p2;
                loadParticle(p2, id_to_index[id2]);
                function(p1, p2, membrane.getDesiredDistance(), membrane.getSpringStrength());
                storeParticle(p1, id_to_index[id1]);
                storeParticle(p2, id_to_index[id2]);
            }
        }

        for(size_t i = 0; i < membrDims[0] - 1; i++){
            for(size_t j = 0; j < membrDims[1]; j++){
                auto id1 = membrane.getMembrNodes()[i][j];
                auto id2 = membrane.getMembrNodes()[i+1][j];

                Particle p1;
                loadParticle(p1, id_to_index[id1]);
                Particle p2;
                loadParticle(p2, id_to_index[id2]);
                function(p1, p2, membrane.getDesiredDistance(), membrane.getSpringStrength());
                storeParticle(p1, id_to_index[id1]);
                storeParticle(p2, id_to_index[id2]);
            }
        }

        for(size_t i = 0; i < membrDims[0]; i++) {
            for (size_t j = 0; j < membrDims[1]; j++) {

                //to top right (thinking about membrane structure)
                if(i+1 < membrDims[0] && j+1 < membrDims[1]){
                    auto id1 = membrane.getMembrNodes()[i][j];
                    auto id2 = membrane.getMembrNodes()[i+1][j+1];

                    Particle p1;
                    loadParticle(p1, id_to_index[id1]);
                    Particle p2;
                    loadParticle(p2, id_to_index[id2]);
                    function(p1, p2, membrane.getDesiredDistance(), membrane.getSpringStrength());
                    storeParticle(p1, id_to_index[id1]);
                    storeParticle(p2, id_to_index[id2]);
                }

                //to bottom right (thinking about membrane structure)
                if(i+1 < membrDims[0] && j-1 >= 0){
                    auto id1 = membrane.getMembrNodes()[i][j];
                    auto id2 = membrane.getMembrNodes()[i+1][j-1];

                    Particle p1;
                    loadParticle(p1, id_to_index[id1]);
                    Particle p2;
                    loadParticle(p2, id_to_index[id2]);
                    function(p1, p2, membrane.getDesiredDistance(), membrane.getSpringStrength());
                    storeParticle(p1, id_to_index[id1]);
                    storeParticle(p2, id_to_index[id2]);
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
                loadParticle(p1, cellItems[i]);
                loadParticle(p2, cellItems[j]);
                function(p1, p2);
                storeParticle(p1, cellItems[i]);
                storeParticle(p2, cellItems[j]);
            }
        }
    }
}

[[maybe_unused]] void ParticleContainer::forAllDistinctCellPairs(
#pragma region param
        void(*fun)(std::vector<double> &force,
                   std::vector<double> &oldForce,
                   std::vector<double> &x,
                   std::vector<double> &v,
                   std::vector<double> &m,
                   std::vector<int> &type,
                   unsigned long count,
                   std::vector<unsigned long> &cell0Items,
                   std::vector<unsigned long> &cell1Items)
#pragma endregion
) {
    io::output::loggers::general->warn(
            "forAllDistinctCellPairs probably wasn't the method you wanted to call you probably wanted to use forAllDistinctCellNeighbours");
    for (unsigned long i = 0; i < cells.size(); i++) {
        for (unsigned long j = i + 1; j < cells.size(); j++) {
            fun(force, oldForce, x, v, m, type, count, cells[i], cells[j]);
        }
    }
}

void ParticleContainer::clearStoreForce() {
    unsigned long size = force.size() / 3;
    for(unsigned long i {0}; i < size; i++) {
        oldForce[3*i + 0] = force[3*i + 0];
        oldForce[3*i + 1] = force[3*i + 1];
        oldForce[3*i + 2] = force[3*i + 2];
        force[3*i + 0] = 0;
        force[3*i + 1] = 0;
        force[3*i + 2] = 0;
    }
}

void ParticleContainer::initAlternativeTaskModel(){
    const unsigned long maxThreads{static_cast<unsigned long>(omp_get_max_threads())};
    alternativeTaskModelCache.clear();
    alternativeTaskModelCache.resize(maxThreads);
    //26 TaskGroups (for the 13 cases*2)
    //every taskGroup has the pairs of CellIndices that are independently doable

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

    /*
    std::array<unsigned int, 3> additionalIncrement{0,0,0};
    for(auto c = 0; c < numCases; c++){ //pun intended
        if(offsets[c][0] == 1){additionalIncrement = std::array<unsigned int, 3>{1,0,0};}
        else if(offsets[c][1] == 1){additionalIncrement = std::array<unsigned int, 3>{0,1,0};}
        else{additionalIncrement = std::array<unsigned int, 3>{0,0,1}*/

    for(auto c = 0; c < numCases; c++){ //pun intended

        #ifdef TASK_ROUND_ROBIN
        constexpr unsigned long roundRobinMolUpdateThreshold = 1'000'000;
        size_t roundRobinAccumulator{0};
        #else
        std::vector<size_t> interactions;
        interactions.resize(maxThreads);
        #endif
        size_t nextIndex{0};

        for(unsigned int x0 = lowerBounds[c][0]; x0 < upperBounds[c][0]; x0++){
            for(unsigned int x1 = lowerBounds[c][1]; x1 < upperBounds[c][1]; x1++){
                for(unsigned int x2 = lowerBounds[c][2]; x2 < upperBounds[c][2]; x2++){
                    auto cell1 = cellIndexFromCellCoordinatesFast(x0, x1, x2);
                    auto cell2 = cellIndexFromCellCoordinatesFast(x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2]);
                    alternativeTaskModelCache[nextIndex].emplace_back(cell1,cell2);
                    SPDLOG_TRACE("Added CellInteraction (({} {} {}), ({} {} {})) to taskBlock {} ", x0, x1, x2, x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2], 2*c+0);

                    #ifdef TASK_ROUND_ROBIN
                    roundRobinAccumulator += cells[cell1].size() * cells[cell2].size();
                    if(roundRobinAccumulator >= roundRobinMolUpdateThreshold){
                        nextIndex = (nextIndex+1)%maxThreads;
                        roundRobinAccumulator = 0;
                    }
                    #else
                    interactions[nextIndex] += cells[cell1].size()*cells[cell2].size();
                    nextIndex = 0;
                    size_t last_count = interactions[nextIndex];
                    for(size_t i = 1; i < maxThreads; i++){
                        if(interactions[i] <= last_count){
                            last_count = interactions[i];
                            nextIndex = i;
                        }
                    }
                    #endif
                }
            }
        }

    }
}

void ParticleContainer::initTaskModel() {
    taskModelCache.clear();
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

    std::array<unsigned int, 3> additionalIncrement{0,0,0};
    for(auto c = 0; c < numCases; c++){ //pun intended
        if(offsets[c][0] == 1){additionalIncrement = std::array<unsigned int, 3>{1,0,0};}
        else if(offsets[c][1] == 1){additionalIncrement = std::array<unsigned int, 3>{0,1,0};}
        else{additionalIncrement = std::array<unsigned int, 3>{0,0,1};}


        //const std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>>& generateDistinctCellNeighbours()
        const unsigned long maxThreads{static_cast<unsigned long>(omp_get_max_threads())};

        #ifdef TASK_ROUND_ROBIN
        constexpr unsigned long roundRobinMolUpdateThreshold = 1'000'000;
        size_t roundRobinAccumulator{0};
        #else
        std::vector<size_t> interactions;
        interactions.resize(maxThreads);
        #endif
        size_t nextIndex{0};

        std::vector<std::vector<std::pair<unsigned long, unsigned long>>> independentTasksBlock{maxThreads, std::vector<std::pair<unsigned long, unsigned long>>{}};

        for(unsigned int x0 = lowerBounds[c][0]; x0 < upperBounds[c][0]; x0+= 1 + additionalIncrement[0]){
            for(unsigned int x1 = lowerBounds[c][1]; x1 < upperBounds[c][1]; x1+= 1 + additionalIncrement[1]){
                for(unsigned int x2 = lowerBounds[c][2]; x2 < upperBounds[c][2]; x2+= 1 + additionalIncrement[2]){

                    auto cell1 = cellIndexFromCellCoordinatesFast(x0, x1, x2);
                    auto cell2 = cellIndexFromCellCoordinatesFast(x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2]);
                    independentTasksBlock[nextIndex].emplace_back(cell1,cell2);
                    SPDLOG_TRACE("Added CellInteraction (({} {} {}), ({} {} {})) to taskBlock {} and job {} ", x0, x1, x2, x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2], 2*c+0, roundRobinIndex);
                    //std::cout<< "Added CellInteraction (("<< x0<<" "<<x1<<" "<< x2 << ") ("<< x0 + offsets[c][0] << " " << x1 + offsets[c][1] <<" "<< x2 + offsets[c][2]<< ")) to taskBlock " << 2*c+0<< " and job " << roundRobinIndex << std::endl;

                    #ifdef TASK_ROUND_ROBIN
                    roundRobinAccumulator += cells[cell1].size() * cells[cell2].size();

                    if(roundRobinAccumulator >= roundRobinMolUpdateThreshold){
                        nextIndex = (nextIndex+1)%maxThreads;
                        roundRobinAccumulator = 0;
                    }
                    #else
                    nextIndex = 0;
                    size_t last_count = interactions[nextIndex];
                    for(size_t i = 1; i < maxThreads; i++){
                        if(interactions[i] <= last_count){
                            last_count = interactions[i];
                            nextIndex = i;
                        }
                    }
                    #endif

                }
            }
        }

        taskModelCache.emplace_back(independentTasksBlock);

        //yes you could cut this down to 2 lines with another helper array but this more verbose version seems much easier to understand
        std::vector<std::vector<std::pair<unsigned long, unsigned long>>> independentTasksBlock2{maxThreads, std::vector<std::pair<unsigned long, unsigned long>>{}};

        #ifdef TASK_ROUND_ROBIN
        nextIndex = 0;
        roundRobinAccumulator = 0;
        #else
        interactions.clear();
        interactions.resize(maxThreads);
        #endif
        for(unsigned int x0 = lowerBounds[c][0] + additionalIncrement[0] ; x0 < upperBounds[c][0]; x0+= 1 + additionalIncrement[0]){
            for(unsigned int x1 = lowerBounds[c][1] + additionalIncrement[1]; x1 < upperBounds[c][1]; x1+= 1 + additionalIncrement[1]){
                for(unsigned int x2 = lowerBounds[c][2] + additionalIncrement[2]; x2 < upperBounds[c][2]; x2+= 1 + additionalIncrement[2]){
                    auto cell1 = cellIndexFromCellCoordinatesFast(x0, x1, x2);
                    auto cell2 = cellIndexFromCellCoordinatesFast(x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2]);
                    independentTasksBlock2[nextIndex].emplace_back(cell1,cell2);
                    SPDLOG_TRACE("Added CellInteraction (({} {} {}), ({} {} {})) to taskBlock {} and job {}", x0, x1, x2, x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2], 2*c+1, roundRobinIndex);
                    //std::cout<< "Added CellInteraction (("<< x0<<" "<<x1<<" "<< x2 << ") ("<< x0 + offsets[c][0] << " " << x1 + offsets[c][1] <<" "<< x2 + offsets[c][2]<< ")) to taskBlock " << 2*c+1<< " and job " << roundRobinIndex << std::endl;

                    #ifdef TASK_ROUND_ROBIN
                    roundRobinAccumulator += cells[cell1].size() * cells[cell2].size();
                    if(roundRobinAccumulator >= roundRobinMolUpdateThreshold){
                        nextIndex = (nextIndex+1)%maxThreads;
                        roundRobinAccumulator = 0;
                        roundRobinAccumulator = 0;
                    }
                    #else
                    nextIndex = 0;
                    size_t last_count = interactions[nextIndex];
                    for(size_t i = 1; i < maxThreads; i++){
                        if(interactions[i] <= last_count){
                            last_count = interactions[i];
                            nextIndex = i;
                        }
                    }
                    #endif

                }
            }
        }
        taskModelCache.emplace_back(independentTasksBlock2);

    }
}



/*
void ParticleContainer::initTaskModel() {
    std::vector<std::vector<std::pair<unsigned long, unsigned long>>> task_group_buffer;
    std::vector<std::pair<unsigned long, unsigned long>> task_buffer;

    //Straight lines ----------------------------------------
    //all pairs in x_0 direction:
    {
        for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
            for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0 + 1, x_1, x_2));
                }
                task_group_buffer.emplace_back(task_buffer);
                task_buffer.clear();
            }
        }
        taskModelCache.emplace_back(task_group_buffer);
        task_group_buffer.clear();

        //all pairs in x_1 direction:
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
            for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0, x_1 + 1, x_2));
                }
                task_group_buffer.emplace_back(task_buffer);
                task_buffer.clear();
            }
        }
        taskModelCache.emplace_back(task_group_buffer);
        task_group_buffer.clear();

        //all pairs in x_2 direction:
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2] - 1; x_2++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0, x_1, x_2 + 1));
                }
                task_group_buffer.emplace_back(task_buffer);
                task_buffer.clear();
            }
        }
        taskModelCache.emplace_back(task_group_buffer);
        task_group_buffer.clear();
    }
    //End of straight lines ---------------------------------------------------

    //"2d-diagonals"------------------------------------------------
    {
        //diagonals lying in the x_0-x_1 plane
        for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
            //diagonals from bottom left to top right
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
                for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0 + 1, x_1 + 1, x_2));
                }
            }
            task_group_buffer.emplace_back(task_buffer);
            task_buffer.clear();
        }
        taskModelCache.emplace_back(task_group_buffer);
        task_group_buffer.clear();

        for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
            //diagonals from top left to bottom right
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
                for (unsigned int x_1 = 1; x_1 < gridDimensions[1]; x_1++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0 + 1, x_1 - 1, x_2));
                }
            }
            task_group_buffer.emplace_back(task_buffer);
            task_buffer.clear();
        }
        taskModelCache.emplace_back(task_group_buffer);
        task_group_buffer.clear();

        //diagonals lying in the x_0-x_2 plane
        for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
            //diagonals from bottom left to top right
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2] - 1; x_2++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0 + 1, x_1, x_2 + 1));
                }
            }
            task_group_buffer.emplace_back(task_buffer);
            task_buffer.clear();
        }
        taskModelCache.emplace_back(task_group_buffer);
        task_group_buffer.clear();
        for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
            //diagonals from top left to bottom right
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
                for (unsigned int x_2 = 1; x_2 < gridDimensions[2]; x_2++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0 + 1, x_1, x_2 - 1));
                }
            }
            task_group_buffer.emplace_back(task_buffer);
            task_buffer.clear();
        }
        taskModelCache.emplace_back(task_group_buffer);
        task_group_buffer.clear();

        //diagonals lying in the x_1-x_2 plane
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
            //diagonals from bottom left to top right
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2] - 1; x_2++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0, x_1 + 1, x_2 + 1));
                }
            }
            task_group_buffer.emplace_back(task_buffer);
            task_buffer.clear();
        }
        taskModelCache.emplace_back(task_group_buffer);
        task_group_buffer.clear();
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
            //diagonals from top left to bottom right
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {
                for (unsigned int x_2 = 1; x_2 < gridDimensions[2]; x_2++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0, x_1 + 1, x_2 - 1));
                }
            }
            task_group_buffer.emplace_back(task_buffer);
            task_buffer.clear();
        }
        taskModelCache.emplace_back(task_group_buffer);
        task_group_buffer.clear();
    }
    //End of "2d diagonals"-----------------------------------------------

    //Start of "3d diagonals"----------------
    {//from bottom front left top back right
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2] - 1; x_2++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0 + 1, x_1 + 1, x_2 + 1));
                }
            }
            task_group_buffer.emplace_back(task_buffer);
            task_buffer.clear();
            taskModelCache.emplace_back(task_group_buffer);
            task_group_buffer.clear();
        }
        //from top front left to bottom back right
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
            for (unsigned int x_1 = 1; x_1 < gridDimensions[1]; x_1++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2] - 1; x_2++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0 + 1, x_1 - 1, x_2 + 1));
                }
            }
            task_group_buffer.emplace_back(task_buffer);
            task_buffer.clear();
            taskModelCache.emplace_back(task_group_buffer);
            task_group_buffer.clear();
        }
        //from bottom back left to top front right
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {
                for (unsigned int x_2 = 1; x_2 < gridDimensions[2]; x_2++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0 + 1, x_1 + 1, x_2 - 1));
                }
            }
            task_group_buffer.emplace_back(task_buffer);
            task_buffer.clear();
            taskModelCache.emplace_back(task_group_buffer);
            task_group_buffer.clear();
        }
        //from top back left to bottom front right
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
            for (unsigned int x_1 = 1; x_1 < gridDimensions[1]; x_1++) {
                for (unsigned int x_2 = 1; x_2 < gridDimensions[2]; x_2++) {
                    task_buffer.emplace_back(cellIndexFromCellCoordinatesFast(x_0, x_1, x_2),
                                             cellIndexFromCellCoordinatesFast(x_0 + 1, x_1 - 1, x_2 - 1));
                }
            }
            task_group_buffer.emplace_back(task_buffer);
            task_buffer.clear();
            taskModelCache.emplace_back(task_group_buffer);
            task_group_buffer.clear();
        }
    }
    //End of "3d diagonals" -----------------
}*/

const std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>>& ParticleContainer::generateDistinctCellNeighbours() {
    initTaskModel();
    return taskModelCache;
}

const std::vector<std::vector<std::pair<unsigned long, unsigned long>>>& ParticleContainer::generateDistinctAlternativeCellNeighbours(){
    initAlternativeTaskModel();
    return alternativeTaskModelCache;
}

#pragma endregion

