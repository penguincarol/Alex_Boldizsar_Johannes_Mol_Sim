//
// Created by johnny on 02.12.22.
//
#include "Particle.h"

#include <vector>

#ifndef PSEMOLDYN_GROUPF_PARTICLECONTAINERREAL_H
#define PSEMOLDYN_GROUPF_PARTICLECONTAINERREAL_H


class ParticleContainerReal {
private:
    double root6_of_2;
    std::vector<Particle> particles; //stores the particles in the particleContainer
    std::vector<unsigned long> activeParticles;
    std::vector<std::vector<unsigned long>> cells;
    std::array<unsigned int, 3> gridDimensions; //stores the number of cells in x- y- and z- direction
    std::array<double, 3> domainSize;
    double x_2_max;
    double x_1_max;
    double x_0_max;
    double x_2_min;
    double x_1_min;
    double x_0_min;
    double r_cutoff;

public:
    /**
     * @brief Construct a new ParticleContainerReal object with no particles stored
     *
     */
    ParticleContainerReal();

    /**
     * @brief Construct a new ParticleContainerReal that stores the given particles
     *
     * @param buffer for all particles. will be added to local storage.
     */
    explicit ParticleContainerReal(const std::vector<Particle> &buffer);

    /**
 * @brief Constructor of ParticleContainerReal that also initializes the cell-structure
 *
 * @param buffer
 * @param domainSize
 * @param r_cutoff
 */
    ParticleContainerReal(const std::vector<Particle>& buffer, std::array<double, 3> domainSize, double r_cutoff);

    /**
     * @brief Constructor of ParticleContainer that also initializes a seemingly two dimensional cell-structure
     *
     * @param buffer
     * @param domainSize
     * @param r_cutoff
     */
    ParticleContainerReal(const std::vector<Particle>& buffer, std::array<double, 2> domainSize, double r_cutoff);

    /**
     * @brief returns the index of the cell in cells corresponding to the coordinates given
     * Example: cellIndexFromCellCoordinates({0,0,0})->0
     * because the cell at position {0,0,0} is stored at index 0 in cells
     * @param coords
     * @return int
     */
    unsigned int cellIndexFromCellCoordinates(std::array<unsigned int, 3> coords);

    /**
     * @brief return the amount of particles stored
     *
     * @return int
     */
    unsigned long size();

    /**
     * Makes sure that every Particle (or every index corresponding to the Particle) is in the
s    * right corresponding cell-vector
     */
    void updateCells();

    /**
     * Removes all particles.
     * */
    void clear();

    /**
     * Get a copy of particle at position @param i
     * */
    Particle getParticle(unsigned long i);
};


#endif //PSEMOLDYN_GROUPF_PARTICLECONTAINERREAL_H
