//
// Created by johnny on 02.12.22.
//
#include "ParticleContainerReal.h"
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

void ParticleContainerReal::runOnData(
        void (*function)(std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &,
                         std::vector<double> &, std::vector<int> &, unsigned long)) {

}
