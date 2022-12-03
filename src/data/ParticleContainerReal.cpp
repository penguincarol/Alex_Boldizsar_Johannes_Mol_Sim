//
// Created by johnny on 02.12.22.
//
#include "ParticleContainerReal.h"
#include "io/output/Logging.h"

#include <Eigen>
#include <numeric>

ParticleContainerReal::ParticleContainerReal(){}

ParticleContainerReal::ParticleContainerReal(const std::vector<Particle> &buffer) {
    particles = std::vector<Particle>(buffer);
}

ParticleContainerReal::ParticleContainerReal(const std::vector<Particle> &buffer, std::array<double, 3> domSize,
                                             double r_cutoff) {
    domainSize = domSize;
    x_2_max = domainSize[2];
    x_1_max = domainSize[1];
    x_0_max = domainSize[0];
    x_2_min = 0;
    x_1_min = 0;
    x_0_min = 0;
    //very explicit casts to absolutely make sure, that the rounding is done correctly
    //this implementation uses shorter grids on the side if the numbers are nasty btw
    std::array<double, 3> helperGridDimensions{std::ceil(domSize[0] / r_cutoff), std::ceil(domSize[1] / r_cutoff),
                                               std::ceil(domSize[2] / r_cutoff)};
    gridDimensions = {(unsigned int) helperGridDimensions[0], (unsigned int) helperGridDimensions[1],
                      (unsigned int) helperGridDimensions[2]};
    //Switch to a different coord-system, where (0,0,0) is the bottom left front corner
    //If we want to use the old coord-system these lines need to get removed and a helper-function should be needed to make the conversion in update-cells
    //and to compute the right array index in this initialization.
    const Eigen::Vector3d offsetCoordConversion{domainSize[0]/2, domainSize[1]/2, domainSize[2]/2};
    for(auto& p: particles){
        p.add_to_X(offsetCoordConversion);
    }

    //i have no idea why i need helper, it should work without it but the compiler doesn't like it
    std::vector<std::vector<unsigned long>> helper(gridDimensions[0] * gridDimensions[1] * gridDimensions[2],
                                                   std::vector<unsigned long>(1)); // TODO fix this
    cells = helper;
    this->r_cutoff = (double) r_cutoff;

    //define which particles are still part of the simulation
    activeParticles.resize(particles.size());
    std::iota(activeParticles.begin(), activeParticles.end(), 0);

    updateCells();

    //halo value
    root6_of_2 = std::pow(std::sqrt(2),3);

}

ParticleContainerReal::ParticleContainerReal(const std::vector<Particle> &buffer, std::array<double, 2> domainSize,
                                 double r_cutoff) :
    ParticleContainerReal::ParticleContainerReal(buffer, {domainSize[0], domainSize[1], r_cutoff}, r_cutoff) {};

unsigned long ParticleContainerReal::size() {
    return particles.size();
}

void ParticleContainerReal::clear() {
    particles.clear();
    activeParticles.clear();
    cells.clear();
}

Particle ParticleContainerReal::getParticle(unsigned long i) {
    return particles[i];
}

std::array<unsigned int, 3> ParticleContainerReal::getGridDimensions() {
    return gridDimensions;
}

void ParticleContainerReal::forAllParticles(const std::function<void(Particle &)> &function) {
    std::for_each(particles.begin(), particles.end(), function);
}
void ParticleContainerReal::forAllParticles(void(*function)(Particle &)) {
    std::for_each(particles.begin(), particles.end(), function);
}

void ParticleContainerReal::forAllPairs(const std::function<void(Particle &, Particle &)> &function) {
    for (unsigned long i{0}; i < particles.size(); i++) {
        for (unsigned long j{i + 1}; j < particles.size(); j++) {
            function(particles[i], particles[j]);
        }
    }
}

unsigned int ParticleContainerReal::cellIndexFromCellCoordinates(std::array<unsigned int, 3> coords) {
//If we decide to change the order of the cells (e.g. by using some fancy 3d space filling curve)
// this method obviously has to be rewritten

//the current version:
//cells[0] corresponds to "cells[0][0][0]"
//cells[1] = "cells[1][0][0]"
//cells[domainSize[0]/r_cutoff] = "cells[0][1][0]"
//cells[(domainSize[0]/r_cutoff)/(domainSize[1]/r_cutoff)] = "cells[0][0][1]"
    unsigned int x0 = std::min(coords[0], gridDimensions[0]-1);
    unsigned int x1 = std::min(coords[1], gridDimensions[1]-1);
    unsigned int x2 = std::min(coords[2], gridDimensions[2]-1);

    return (x0 +
            x1 * gridDimensions[0] +
            x2 * gridDimensions[0] * gridDimensions[1]
    );
}

void ParticleContainerReal::forAllPairsInSameCell(const std::function<void(Particle &p1, Particle &p2)>& function){
    for (std::vector<unsigned long> &cellItems : cells){
        for(unsigned long i=0; i<cellItems.size(); i++){
            for(unsigned long j=i+1;j<cellItems.size(); j++){
                Particle p1, p2;//TODO 0 check this
                function(particles[cellItems[i]], particles[cellItems[j]]);
            }
        }
    }
}

void ParticleContainerReal::forAllPairsInNeighbouringCell(const std::function<void(Particle &p1, Particle &p2)>& function){

    //Implementation2:
    //basically every code snippet occurs three times right here because every dimension needs to be the "free variable" for every case once
    //but actually this seems more robust than making some fancy "iterate over all possible variable distribution"-thingies

    //Straight lines ----------------------------------------
    //all pairs in x_0 direction:
    for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
        for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0 + 1, x_1, x_2})]){
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0 + 1, x_1, x_2);
            }
        }
    }
    //all pairs in x_1 direction:
    for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
        for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0, x_1 + 1, x_2})]){
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0, x_1 + 1, x_2);
            }
        }
    }
    //all pairs in x_2 direction:
    for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
        for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
            for (unsigned int x_2 = 0; x_2 < gridDimensions[2] - 1; x_2++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2 + 1})]){
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0, x_1, x_2 + 1);
            }
        }
    }
    //End of straight lines ---------------------------------------------------

    //"2d-diagonals"------------------------------------------------
    //diagonals lying in the x_0-x_1 plane
    for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
        //diagonals from bottom left to top right
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0 + 1, x_1 + 1,x_2})]){ //check with the neighbour that is one to the right and one above you
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0 + 1, x_1 + 1,x_2);
            }
        }
        //diagonals from top left to bottom right
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
            for (unsigned int x_1 = 1; x_1 < gridDimensions[1]; x_1++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0 + 1, x_1 - 1,x_2})]){ //(check with the neighbour that is one to the right and one below you)
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0 + 1, x_1 - 1,x_2);
            }
        }
    }
    //diagonals lying in the x_0-x_2 plane
    for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
        //diagonals from bottom left to top right
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
            for (unsigned int x_2 = 0; x_2 < gridDimensions[2] - 1; x_2++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0 + 1, x_1, x_2 +1})]){ //check with the neighbour that is one to the right and one above you
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0 + 1, x_1, x_2 +1);
            }
        }
        //diagonals from top left to bottom right
        for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
            for (unsigned int x_2 = 1; x_2 < gridDimensions[2]; x_2++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0 + 1, x_1, x_2 -
                                                                                                 1})]){ //(check with the neighbour that is one to the right and one below you)
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0 + 1, x_1, x_2 -1);
            }
        }
    }
    //diagonals lying in the x_1-x_2 plane
    for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
        //diagonals from bottom left to top right
        for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {
            for (unsigned int x_2 = 0; x_2 < gridDimensions[2] - 1; x_2++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0, x_1 + 1, x_2 +
                                                                                                 1})]){ //(check with the neighbour that is one to the right and one below you)
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0, x_1 + 1, x_2 + 1);
            }
        }
        //diagonals from top left to bottom right
        for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {
            for (unsigned int x_2 = 1; x_2 < gridDimensions[2]; x_2++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0, x_1 + 1, x_2 -
                                                                                                 1})]){ //(check with the neighbour that is one to the right and one below you)
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0, x_1 + 1, x_2 -1);
            }
        }
    }
    //End of "2d diagonals"-----------------------------------------------

    //Start of "3d diagonals"----------------
    //from bottom front left top back right
    for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
        for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {
            for (unsigned int x_2 = 0; x_2 < gridDimensions[2] - 1; x_2++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0 + 1, x_1 + 1, x_2 + 1})]){
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0 + 1, x_1 + 1, x_2 + 1);
                //std::cout<<"(" << x_0 << ", " << x_1 << ", " << x_2 << ") interacted with (" << x_0+1 << ", " << x_1+1 << ", " << x_2+1 << ")\n";
            }
        }
    }
    //from top front left to bottom back right
    for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
        for (unsigned int x_1 = 1; x_1 < gridDimensions[1]; x_1++) {
            for (unsigned int x_2 = 0; x_2 < gridDimensions[2] - 1; x_2++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0 + 1, x_1 - 1, x_2 + 1})]){
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0 + 1, x_1 - 1, x_2 + 1);
            }
        }
    }
    //from bottom back left to top front right
    for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
        for (unsigned int x_1 = 0; x_1 < gridDimensions[1] - 1; x_1++) {
            for (unsigned int x_2 = 1; x_2 < gridDimensions[2]; x_2++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0 + 1, x_1 + 1, x_2 - 1})]){
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0 + 1, x_1 + 1, x_2 - 1);
            }
        }
    }
    //from top back left to bottom front right
    for (unsigned int x_0 = 0; x_0 < gridDimensions[0] - 1; x_0++) {
        for (unsigned int x_1 = 1; x_1 < gridDimensions[1]; x_1++) {
            for (unsigned int x_2 = 1; x_2 < gridDimensions[2]; x_2++) {

                for(unsigned long index0 : cells[cellIndexFromCellCoordinates({x_0, x_1, x_2})]){
                    for(unsigned long index1 : cells[cellIndexFromCellCoordinates({x_0 + 1, x_1 - 1, x_2 - 1})]){
                        function(particles[index0], particles[index1]);
                    }
                }
                io::output::loggers::general->trace("Cell ({} {} {}) interacted with ({} {} {})", x_0, x_1, x_2,
                                                    x_0 + 1, x_1 - 1, x_2 - 1);
            }
        }
    }
}

void ParticleContainerReal::updateCells() {
    //I am doing an implementation that works first and then i figure out if there is a better way
    //than deciding for every particle in every iteration once again

    //by the way: is there a way to advice a vector not to shrink? i can't find it with like.. 10 mins of googling
    io::output::loggers::general->trace("updateCells called");
    for (auto& cell: cells) {
        cell.clear();
    }
    for (unsigned int i: activeParticles) {
        //i am intentionally rounding down with casts from double to unsigned int
        std::array<unsigned int, 3> cellCoordinate{(unsigned int) (particles[i].getX()[0] / r_cutoff),
                                                   (unsigned int) (particles[i].getX()[1] / r_cutoff),
                                                   (unsigned int) (particles[i].getX()[2]  / r_cutoff)};
        this->cells[cellIndexFromCellCoordinates(cellCoordinate)].emplace_back(i);
    }
}

void ParticleContainerReal::forAllCells(void (*fun)(std::vector<double> &force,
                             std::vector<double> &oldForce,
                             std::vector<double> &x,
                             std::vector<double> &v,
                             std::vector<double> &m,
                             std::vector<int> &type,
                             unsigned long count,
                             std::vector<unsigned long>& cellItems)){
    {
        for (auto &cellItems: cells) {
            std::vector<double> force{};
            std::vector<double> oldForce{};
            std::vector<double> x{};
            std::vector<double> v{};
            std::vector<double> m{};
            std::vector<int> type{};
            unsigned long count{0};
            fillFlattenedVectors(particles, force, oldForce, x, v, m, type, count);

            fun(force, oldForce, x, v, m, type, count, cellItems);

            writeBackFlattenedVectors(force, oldForce, x, v, m, type, count);
        }
    }
}

void ParticleContainerReal::forAllDistinctCellPairs(void (*fun)(std::vector<double> &force,
                                         std::vector<double> &oldForce,
                                         std::vector<double> &x,
                                         std::vector<double> &v,
                                         std::vector<double> &m,
                                         std::vector<int> &type,
                                         unsigned long count,
                                         std::vector<unsigned long>& cell0Items,
                                         std::vector<unsigned long>& cell1Items)){
    io::output::loggers::general->warn(
            "forAllDistinctCellPairs probably wasn't the method you wanted to call you probably wanted to use forAllDistinctCellNeighbours");
    for (unsigned long i = 0; i < cells.size(); i++) {
        for (unsigned long j = i + 1; j < cells.size(); j++) {
            std::vector<double> force{};
            std::vector<double> oldForce{};
            std::vector<double> x{};
            std::vector<double> v{};
            std::vector<double> m{};
            std::vector<int> type{};
            unsigned long count{0};
            fillFlattenedVectors(particles, force, oldForce, x, v, m, type, count);

            fun(force, oldForce, x, v, m, type, count, cells[i], cells[j]);

            writeBackFlattenedVectors(force, oldForce, x, v, m, type, count);
        }
    }
}

void ParticleContainerReal::fillFlattenedVectors(const std::vector<Particle> &buffer,
                                                 std::vector<double> &force,
                                                 std::vector<double> &oldForce,
                                                 std::vector<double> &x,
                                                 std::vector<double> &v,
                                                 std::vector<double> &m,
                                                 std::vector<int> &type,
                                                 unsigned long &count) {
    count = buffer.size();
    force.resize(count * 3);
    oldForce.resize(count * 3);
    x.resize(count * 3);
    v.resize(count * 3);
    m.resize(count);
    type.resize(count);

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
    }
}

void ParticleContainerReal::writeBackFlattenedVectors(std::vector<double> &force,
                               std::vector<double> &oldForce,
                               std::vector<double> &x,
                               std::vector<double> &v,
                               std::vector<double> &m,
                               std::vector<int> &type,
                               unsigned long &count){
    for(unsigned int i{0}; i < count; i++){
        Particle& p{particles[i]};
        p.setX({x[3*i], x[3*i +1], x[3*i +2]});
        p.setF({force[3*i], force[3*i +1], force[3*i +2]});
        p.setOldF({oldForce[3*i], oldForce[3*i +1], oldForce[3*i +2]});
        p.setV({v[3*i], v[3*i +1], v[3*i +2]});
        p.setM(m[i]);
        p.setType(type[i]);
    }
}

void ParticleContainerReal::runOnData(
        void (*function)(std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &,
                         std::vector<double> &, std::vector<int> &, unsigned long)) {
    std::vector<double> force{};
    std::vector<double> oldForce{};
    std::vector<double> x{};
    std::vector<double> v{};
    std::vector<double> m{};
    std::vector<int> type{};
    unsigned long count{0};
    fillFlattenedVectors(particles, force, oldForce, x, v, m, type, count);
    function(force, oldForce, x, v, m, type, count);
    writeBackFlattenedVectors(force, oldForce, x, v, m, type, count);
}

void ParticleContainerReal::runOnData(void (*fun)(std::vector<double> &force,
                           std::vector<double> &oldForce,
                           std::vector<double> &x,
                           std::vector<double> &v,
                           std::vector<double> &m,
                           std::vector<int> &type,
                           unsigned long count,
                           std::vector<std::vector<unsigned long>>& cells)){
    std::vector<double> force{};
    std::vector<double> oldForce{};
    std::vector<double> x{};
    std::vector<double> v{};
    std::vector<double> m{};
    std::vector<int> type{};
    unsigned long count{0};
    fillFlattenedVectors(particles, force, oldForce, x, v, m, type, count);
    fun(force, oldForce, x, v, m, type, count, cells);
    writeBackFlattenedVectors(force, oldForce, x, v, m, type, count);
}

void ParticleContainerReal::runOnData(const std::function<void(std::vector<double> &force,
                                        std::vector<double> &oldForce,
                                        std::vector<double> &x,
                                        std::vector<double> &v,
                                        std::vector<double> &m,
                                        std::vector<int> &type,
                                        unsigned long count)>& function){
    std::vector<double> force{};
    std::vector<double> oldForce{};
    std::vector<double> x{};
    std::vector<double> v{};
    std::vector<double> m{};
    std::vector<int> type{};
    unsigned long count{0};
    fillFlattenedVectors(particles, force, oldForce, x, v, m, type, count);
    function(force, oldForce, x, v, m, type, count);
    writeBackFlattenedVectors(force, oldForce, x, v, m, type, count);
}