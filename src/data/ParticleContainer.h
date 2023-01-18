#pragma once

#include "ParticleContainer.h"
#include "Particle.h"
#include "sim/physics/bounds/types.h"
#include "io/output/Logging.h"
#include "Membrane.h"

#include <vector>
#include <array>
#include <string>
#include <iterator>
#include <functional>
#include <unordered_set>

/**
 * @brief wrapper class that stores and manages access to the particles
 *      The given implementation is a wrapper class around a std::vector<Particle>
 *      Another implementation might decide to use a different underlying structure.
 */
class ParticleContainer {
public:
    /**
     * Is in itself a 3d vector, but can provide access to a 3d sub-vector with a certain inset
     * s.t. the inner vector seems as if there is no outer vector.
     * The inset is defined by offsetXI. This offset is applied to either side individually.
     * Resulting in a total difference of 2*offset;
     * */
    class VectorCoordWrapper {
    private:
        std::vector<std::vector<unsigned long>> p_data;
        unsigned long p_dimX0;
        unsigned long p_dimX1;
        unsigned long p_dimX2;
        unsigned long p_offsetX0;
        unsigned long p_offsetX1;
        unsigned long p_offsetX2;
        unsigned long p_eDimX0;
        unsigned long p_eDimX1;
        unsigned long p_eDimX2;

    public:
        struct Iterator {
            using iterator_category = std::forward_iterator_tag;
            using difference_type = std::ptrdiff_t;
            using value_type = std::vector<unsigned long>;
            using pointer = std::vector<unsigned long> *;
            using reference = std::vector<unsigned long> &;

        private:
            unsigned long p_index;
            VectorCoordWrapper &p_vector;
        public:
            explicit Iterator(unsigned long i, VectorCoordWrapper &vector) : p_index(i), p_vector(vector) {}

            reference operator*() const { return p_vector.getInner(p_index); }

            pointer operator->() { return &p_vector.getInner(p_index); }

            // Prefix increment
            Iterator &operator++() {
                p_index++;
                return *this;
            }

            // Postfix increment
            Iterator operator++(int) {
                Iterator tmp = *this;
                ++(*this);
                return tmp;
            }

            friend bool operator==(const Iterator &a, const Iterator &b) { return a.p_index == b.p_index; };

            friend bool operator!=(const Iterator &a, const Iterator &b) { return a.p_index != b.p_index; };
        };

        explicit VectorCoordWrapper(unsigned long dimX0 = 0, unsigned long dimX1 = 0, unsigned long dimX2 = 0,
                                    unsigned long dx0 = 1, unsigned long dx1 = 1, unsigned long dx2 = 1) :
                p_dimX0(dimX0), p_dimX1(dimX1), p_dimX2(dimX2),
                p_offsetX0(dx0), p_offsetX1(dx1), p_offsetX2(dx2),
                p_eDimX0(dimX0 - 2 * dx0), p_eDimX1(dimX1 - 2 * dx1), p_eDimX2(dimX2 - 2 * dx2) {
            p_data.resize(dimX0 * dimX1 * dimX2);
        }

        VectorCoordWrapper(VectorCoordWrapper &o) noexcept: p_dimX0(o.p_dimX0), p_dimX1(o.p_dimX1), p_dimX2(o.p_eDimX2),
                                                            p_offsetX0(o.p_offsetX0), p_offsetX1(o.p_offsetX1),
                                                            p_offsetX2(o.p_offsetX2),
                                                            p_eDimX0(o.p_eDimX0), p_eDimX1(o.p_eDimX1),
                                                            p_eDimX2(o.p_eDimX2) {
            p_data.resize(o.p_data.size());
        }

        /**
         * Returns the vector at the global index.
         * */
        [[maybe_unused]] std::vector<unsigned long> &getOuter(unsigned long index) {
            return p_data[index];
        }

        /**
         * Returns the vector at the given coordinates
         * */
        [[maybe_unused]] std::vector<unsigned long> &getOuter(unsigned long x0, unsigned long x1, unsigned long x2) {
            return p_data[x0 + x1 * p_dimX0 + x2 * p_dimX0 * p_dimX1];
        }

        /**
         * Returns the vector at the given inner coordinates. \n
         * Inner to global coord conversion follows: \n
         * global: g = x0 + x1 * d0 + x2 * d0*d1 \n
         * inner : i = y0 + y1 * e0 + y2 * e0*e1 = \n
         *             (x0 + 1) + (x1 + 1) * e0 + (x2 + 1) * e0*e1 = \n
         *             (x0 + 1) + (x1 + 1) * (d0 - 2*o0) + (x2 + 1) * (d0 - 2*o0)*(d1 - 2*o1) = \n
         *             x0 + x1 * d0 + x2 * d0*d1 - 2x2d0o1 - 2x2d1o0 + 4x2o0o1 + d0d1 - 2d0o1 - 2d1o0 + 4o0o1 + 1 - x1 * 2*o0 + d0 - 2*o0 = \n
         *  ==>   g = i - (- 2x2d0o1 - 2x2d1o0 + 4x2o0o1 + d0d1 - 2d0o1 - 2d1o0 + 4o0o1 + 1 - x1 * 2*o0 + d0 - 2*o0) = \n
         *            i + 2(y2 - 1)d0o1 + 2(y2 - 1)d1o0 - 4(y2 - 1)o0o1 - d0d1 + 2d0o1 + 2d1o0 - 4o0o1 - 1 + 2(y1 - 1)o0 - d0 + 2*o0) = \n
         *            i + 2y2d0o1 - 2d0o1 + 2y2d1o0 - 2d1o0 - 4y2o0o1 + 4o0o1 - d0d1 + 2d0o1 + 2d1o0 - 4o0o1 - 1 + 2y1o0 - 2o0 - d0 + 2o0 = \n
         *            i + 2y2d0o1 + 2y2d1o0 - 4y2o0o1 - d0d1 - 1 + 2y1o0 - d0 = \n
         *  @param y0 x0 coord in inner coord system
         *  @param y1 x1 coord in inner coord system
         *  @param y2 x2 coord in inner coord system
         * */
        [[maybe_unused]] std::vector<unsigned long> &getInner(unsigned long y0, unsigned long y1, unsigned long y2) {
            unsigned long g = (y0 + 1) + (y1 + 1) * p_dimX0 + (y2 + 1) * p_dimX0 * p_dimX1;
            return p_data[g];
//            unsigned long i = y0 + y1 * p_eDimX0 + y2 * p_eDimX0 * p_eDimX1;
//            unsigned long g = i + 2 * y2 * p_dimX0 * p_offsetX1 + 2 * y2 * p_dimX1 * p_offsetX0
//                                - 4 * y2 * p_offsetX0 * p_offsetX1 - p_dimX0 * p_dimX1 - 1
//                                + 2 * y1 * p_offsetX0 - p_dimX0;
//            return p_data[g];
        }

        /**
         * Returns the vector at the given inner coordinates. \n
         * Inner to global coord conversion follows: \n
         * global: g = x0 + x1 * d0 + x2 * d0*d1 \n
         * inner : i = y0 + y1 * e0 + y2 * e0*e1 = \n
         *             (x0 + 1) + (x1 + 1) * e0 + (x2 + 1) * e0*e1 = \n
         *             (x0 + 1) + (x1 + 1) * (d0 - 2*o0) + (x2 + 1) * (d0 - 2*o0)*(d1 - 2*o1) = \n
         *             x0 + x1 * d0 + x2 * d0*d1 - 2x2d0o1 - 2x2d1o0 + 4x2o0o1 + d0d1 - 2d0o1 - 2d1o0 + 4o0o1 + 1 - x1 * 2*o0 + d0 - 2*o0 = \n
         *  ==>   g = i - (- 2x2d0o1 - 2x2d1o0 + 4x2o0o1 + d0d1 - 2d0o1 - 2d1o0 + 4o0o1 + 1 - x1 * 2*o0 + d0 - 2*o0) = \n
         *            i + 2(y2 - 1)d0o1 + 2(y2 - 1)d1o0 - 4(y2 - 1)o0o1 - d0d1 + 2d0o1 + 2d1o0 - 4o0o1 - 1 + 2(y1 - 1)o0 - d0 + 2*o0) = \n
         *            i + 2y2d0o1 - 2d0o1 + 2y2d1o0 - 2d1o0 - 4y2o0o1 + 4o0o1 - d0d1 + 2d0o1 + 2d1o0 - 4o0o1 - 1 + 2y1o0 - 2o0 - d0 + 2o0 = \n
         *            i + 2y2d0o1 + 2y2d1o0 - 4y2o0o1 - d0d1 - 1 + 2y1o0 - d0 = \n
         *  @param y0 x0 coord in inner coord system
         *  @param y1 x1 coord in inner coord system
         *  @param y2 x2 coord in inner coord system
         * */
        [[maybe_unused]] std::vector<unsigned long> &getInner(unsigned long i) {
            unsigned long y0;
            unsigned long y1;
            unsigned long y2;
            y0 = i % p_eDimX0;
            y1 = ((i - y0) % (p_eDimX0 * p_eDimX1)) / p_eDimX0;
            y2 = (i - y0 - y1 * p_eDimX0) / (p_eDimX0 * p_eDimX1);
            return getInner(y0, y1, y2);
//            y1 = (i % (p_eDimX0 * p_eDimX1)) / p_eDimX1;
//            y2 = i / (p_eDimX0 * p_eDimX1);
//            unsigned long g = i + 2 * y2 * p_dimX0 * p_offsetX1 + 2 * y2 * p_dimX1 * p_offsetX0
//                              - 4 * y2 * p_offsetX0 * p_offsetX1 - p_dimX0 * p_dimX1 - 1
//                              + 2 * y1 * p_offsetX0 - p_dimX0;
//            return p_data[g];
        }

        /**
         * delegates to getInner(i)
         * Operates on the inner 3d vector
         * */
        [[maybe_unused]] std::vector<unsigned long> &operator[](unsigned long i) {
            return getInner(i);
        }

        /**
         * Operates only on the inner 3d vector
         * */
        Iterator begin() { return Iterator(0, *this); }

        /**
         * Operates only on the inner 3d vector
         * */
        Iterator end() { return Iterator(p_eDimX0 * p_eDimX1 * p_eDimX2, *this); }

        /**
         * Clears all cells, also halo
         * */
        void clear() {
            p_data.clear();
        }

        /**
         * Clears all outer cells
         * */
        void clearOuter() {
            // l r
            for (unsigned long o0 = 0; o0 < p_offsetX0; o0++) {
                for (unsigned long x1 = 0; x1 < p_dimX1; x1++) {
                    for (unsigned long x2 = 0; x2 < p_dimX2; x2++) {
                        getOuter(0 + o0, x1, x2).clear();
                        getOuter(p_dimX0 - 1 - o0, x1, x2).clear();
                    }
                }
            }
            // b t
            for (unsigned long o1 = 0; o1 < p_offsetX1; o1++) {
                for (unsigned long x0 = 0; x0 < p_dimX0; x0++) {
                    for (unsigned long x2 = 0; x2 < p_dimX2; x2++) {
                        getOuter(x0, 0 + o1, x2).clear();
                        getOuter(x0, p_dimX1 - 1 - o1, x2).clear();
                    }
                }
            }
            // f r
            for (unsigned long o2 = 0; o2 < p_offsetX2; o2++) {
                for (unsigned long x0 = 0; x0 < p_dimX0; x0++) {
                    for (unsigned long x1 = 0; x1 < p_dimX1; x1++) {
                        getOuter(x0, x1, 0 + o2).clear();
                        getOuter(x0, x1, p_dimX2 - 1 - o2).clear();
                    }
                }
            }
        }

        /**
         * Returns the size of the inner 3d vector
         * */
        unsigned long size() {
            return p_eDimX0 * p_eDimX1 * p_eDimX2;
        }

        /**
         * Stores value into cell at given coordinates.
         * Will check if it is already contained in cell.
         * Will not store a duplicate.
         * Coords are in global coord system.
         * */
        void store(unsigned long x0, unsigned long x1, unsigned long x2, unsigned long val) {
            x0 = std::min(x0, p_dimX0 - 1);
            x1 = std::min(x1, p_dimX1 - 1);
            x2 = std::min(x2, p_dimX2 - 1);
            auto &cell = getOuter(x0, x1, x2);
            if (std::find(cell.begin(), cell.end(), val) != cell.end()) return;
            cell.emplace_back(val);
        }

        /**
         * Returns an inner cell using by addressing it with the global coord system.
         * If the specified point is outside of the inner 3d vector or OOB then an empty vector is returned
         * */
        std::vector<unsigned long> &getInnerGlobal(unsigned long x0, unsigned long x1, unsigned long x2) {
            static std::vector<unsigned long> empty;
            if (x0 < p_dimX0 - p_offsetX0 - p_eDimX0 || x0 > p_offsetX0 + p_eDimX0 - 1) return empty;
            if (x1 < p_dimX1 - p_offsetX1 - p_eDimX1 || x1 > p_offsetX1 + p_eDimX1 - 1) return empty;
            if (x2 < p_dimX2 - p_offsetX2 - p_eDimX2 || x2 > p_offsetX2 + p_eDimX2 - 1) return empty;
            return getOuter(x0, x1, x2);
        }

        /**
         * Removes val from cell i
         * @param i inner dense coordinate
         * @param val item to remove
         * */
        void removeAt(unsigned long i, unsigned long val) {
            auto &cell = getInner(i);
            cell.erase(std::remove_if(cell.begin(), cell.end(), [&](const auto item) { return item == val; }),
                       cell.end());
        }
    };

private:
    /**
     * Amount of indices, to be padded with.
     * Resulting padding size is: padding_count * 3 * 8 in Byte.
     * E.G. with padding count 6 => padded by 144 Byte
     * */
    static constexpr unsigned long padding_count = 1;
    double root6_of_2;
    std::vector<Particle> particles;
    unsigned long count;
    /**contains particle IDs*/
    std::vector<unsigned long> activeParticles;
    VectorCoordWrapper cells;
    std::array<unsigned int, 3> gridDimensions; //stores the number of cells in x- y- and z- direction
    std::array<double, 3> domainSize;
    std::vector<Membrane> membranes;  //stores the membranes created membranes as defined in the Membrane-class
    double x_2_max;
    double x_1_max;
    double x_0_max;
    double x_2_min;
    double x_1_min;
    double x_0_min;
    double r_cutoff;
    std::unordered_map<unsigned long, unsigned long> id_to_index;
    std::unordered_map<unsigned long, unsigned long> index_to_id;
    bool eOMP;

    /**
     * Swaps the position of the given particles. Addressed using particle IDs. Location will be looked up in id_to_index map.
     * */
    void swap(unsigned long id0, unsigned long id1);

    /**
     * Moves particle at index indexSrc to index indexDst. id is the ID of the particle to be moved.
     * */
    void move(unsigned long indexSrc, unsigned indexDst, unsigned long id);

public:
    /**
     * @brief Construct a new Particle Container object with no particles stored
     *
     */
    ParticleContainer();

    /**
     * @brief Construct a new Particle Container that stores the given particles
     *
     * @param buffer for all particles. will be added to local storage.
     */
    explicit ParticleContainer(const std::vector<Particle> &buffer);

    /**
     * @brief Constructor of ParticleContainer that also initializes the cell-structure
     *
     * @param buffer
     * @param domainSize
     * @param r_cutoff
     */
    ParticleContainer(const std::vector<Particle> &buffer, std::array<double, 3> domainSize, double r_cutoff, const std::vector<Membrane>& membranes = {}, bool fEOMP = false);

    /**
     * @brief Constructor of ParticleContainer that also initializes a seemingly two dimensional cell-structure
     *
     * @param buffer
     * @param domainSize
     * @param r_cutoff
     */
    ParticleContainer(const std::vector<Particle> &buffer, std::array<double, 2> domainSize, double r_cutoff, const std::vector<Membrane>& membranes = {}, bool fEOMP = false);

    /**
     * @brief returns the index of the cell in cells corresponding to the coordinates given. Performs NO bounds checks!
     * Example: cellIndexFromCellCoordinates({0,0,0})->0
     * because the cell at position {0,0,0} is stored at index 0 in cells
     * @param coords
     * @return int
     */
    unsigned int cellIndexFromCellCoordinatesFast(unsigned int x0, unsigned int x1, unsigned int x2);

    /**
     * @brief returns the index of the cell in cells corresponding to the coordinates given. Performs bounds checks
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
     * @brief Returns the amount of the active particles
     * */
    unsigned long activeSize();

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

    /**
     * @return reference to internally stored membranes
     */
    std::vector<Membrane>& getMembranes();

    /**
     * Moves all forces to the oldForces buffer
     * */
     void clearStoreForce();

    /**
     * @brief Generates all halo particles, that are not stored yet
     * */
    template<sim::physics::bounds::side S>
    void
    forEachParticleHaloPairInSide(const double sigma, const std::function<void(Particle &, Particle &)> &function) {
        double maxBorderDistance = root6_of_2 * sigma / 2;
        if constexpr (S == sim::physics::bounds::side::left) {
            // left x = x_0 = min
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_left = cells[cellIndexFromCellCoordinates({0, x_1, x_2})];
                    for (auto i: cell_indices_left) {
                        Particle& p_real = particles[i];
                        auto& x = p_real.getX();
                        if (x[0] != 0 && x[0] < maxBorderDistance) {
                            Particle p_halo = p_real;
                            p_halo.setX({(-1) * x[0], x[1], x[2]});
                            function(p_real, p_halo);
                        }
                    }
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::right) {
            // right x = x_0 = max and left x = x_0 = min
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_right = cells[cellIndexFromCellCoordinates({gridDimensions[0] - 1, x_1, x_2})];
                    for (auto i: cell_indices_right) {
                        Particle& p_real = particles[i];
                        auto& x = p_real.getX();
                        double distance = domainSize[0] - x[0];
                        if (distance < maxBorderDistance) {
                            Particle p_halo = p_real;
                            p_halo.add_to_X({distance * 2, 0, 0});
                            function(p_real, p_halo);
                        }
                    }
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::bottom) {
            // top y = x_1 = max and bottom y = x_1 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_bot = cells[cellIndexFromCellCoordinates({x_0, 0, x_2})];
                    for (auto i: cell_indices_bot) {
                        Particle& p_real = particles[i];
                        auto& x = p_real.getX();
                        if (x[1] != 0 && x[1] < maxBorderDistance) {
                            Particle p_halo = p_real;
                            p_halo.setX({x[0], (-1) * x[1], x[2]});

                            function(p_real, p_halo);
                        }
                    }
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::top) {
            // top y = x_1 = max and bottom y = x_1 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_top = cells[cellIndexFromCellCoordinates({x_0, gridDimensions[1] - 1, x_2})];
                    for (auto i: cell_indices_top) {
                        Particle& p_real = particles[i];
                        auto& x = p_real.getX();
                        double distance = domainSize[1] - x[1];
                        if (distance < maxBorderDistance) {
                            Particle p_halo = p_real;
                            p_halo.add_to_X({0, distance * 2, 0});

                            function(p_real, p_halo);
                        }
                    }
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::front) {
            // back z = x_2 = max and front z = x_2 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                    auto &cell_indices_front = cells[cellIndexFromCellCoordinates({x_0, x_1, 0})];
                    for (auto i: cell_indices_front) {
                        Particle& p_real = particles[i];
                        auto& x = p_real.getX();
                        if (x[2] != 0 && x[2] < maxBorderDistance) {
                            Particle p_halo = p_real;
                            p_halo.setX({x[0], x[1], (-1) * x[2]});

                            function(p_real, p_halo);
                        }
                    }
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::rear) {
            // back z = x_2 = max and front z = x_2 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                    auto &cell_indices_back = cells[cellIndexFromCellCoordinates({x_0, x_1, gridDimensions[2] - 1})];
                    for (auto i: cell_indices_back) {
                        Particle& p_real = particles[i];
                        auto& x = p_real.getX();
                        double distance = domainSize[2] - x[2];
                        if (distance < maxBorderDistance) {
                            Particle p_halo = p_real;
                            p_halo.add_to_X({0, 0, distance * 2});

                            function(p_real, p_halo);
                        }
                    }
                }
            }
        }
    }

    /**
     * Pops the indices of all particles, that are outside of the domain.
     * Will only use the correct branch. This is determined at compile time.
     * Indices of particles outside the domain are removed of the cells data structure and stored into output.
     * */
    template<sim::physics::bounds::side S>
    void popExternalParticles(std::unordered_set<unsigned long> &output) {
        if constexpr (S == sim::physics::bounds::side::left) {
            // right x = x_0 = max and left x = x_0 = min
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_left = cells[cellIndexFromCellCoordinates({0, x_1, x_2})];
                    cell_indices_left.erase(
                            std::remove_if(cell_indices_left.begin(), cell_indices_left.end(), [&](auto i) {
                                if (particles[i].getX()[0] < x_0_min) {
                                    output.emplace(i);
                                    return true;
                                }
                                return false;
                            }),
                            cell_indices_left.end()
                    );
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::right) {
            // right x = x_0 = max and left x = x_0 = min
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_right = cells[cellIndexFromCellCoordinates({gridDimensions[0] - 1, x_1, x_2})];
                    cell_indices_right.erase(
                            std::remove_if(cell_indices_right.begin(), cell_indices_right.end(), [&](auto i) {
                                if (particles[i].getX()[0] > x_0_max) {
                                    output.emplace(i);
                                    return true;
                                }
                                return false;
                            }),
                            cell_indices_right.end()
                    );
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::bottom) {
            // top y = x_1 = max and bottom y = x_1 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_bot = cells[cellIndexFromCellCoordinates({x_0, 0, x_2})];
                    cell_indices_bot.erase(
                            std::remove_if(cell_indices_bot.begin(), cell_indices_bot.end(), [&](auto i) {
                                if (particles[i].getX()[1] < x_1_min) {
                                    output.emplace(i);
                                    return true;
                                }
                                return false;
                            }),
                            cell_indices_bot.end()
                    );
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::top) {
            // top y = x_1 = max and bottom y = x_1 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_top = cells[cellIndexFromCellCoordinates({x_0, gridDimensions[1] - 1, x_2})];
                    cell_indices_top.erase(
                            std::remove_if(cell_indices_top.begin(), cell_indices_top.end(), [&](auto i) {
                                if (particles[i].getX()[1] > x_1_max) {
                                    output.emplace(i);
                                    return true;
                                }
                                return false;
                            }),
                            cell_indices_top.end()
                    );
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::front) {
            // back z = x_2 = max and front z = x_2 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                    auto &cell_indices_front = cells[cellIndexFromCellCoordinates({x_0, x_1, 0})];
                    cell_indices_front.erase(
                            std::remove_if(cell_indices_front.begin(), cell_indices_front.end(), [&](auto i) {
                                if (particles[i].getX()[2] < x_2_min) {
                                    output.emplace(i);
                                    return true;
                                }
                                return false;
                            }),
                            cell_indices_front.end()
                    );
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::rear) {
            // back z = x_2 = max and front z = x_2 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                    auto &cell_indices_back = cells[cellIndexFromCellCoordinates({x_0, x_1, gridDimensions[2] - 1})];
                    cell_indices_back.erase(
                            std::remove_if(cell_indices_back.begin(), cell_indices_back.end(), [&](auto i) {
                                if (particles[i].getX()[2] > x_2_max) {
                                    output.emplace(i);
                                    return true;
                                }
                                return false;
                            }),
                            cell_indices_back.end()
                    );
                }
            }
        }
    }

    /**
     * Removes the specified indices from the active list. For this to have effect updateCells has to be called.
     * */
    void deactivateParticles(std::unordered_set<unsigned long> &indices);

    /**
     * moves the particles at the given indices from one bound side to the opposing one.
     * will move all cell indices as well
     * */
    template<sim::physics::bounds::side S>
    void moveExternalParticles(std::unordered_set<unsigned long> &indices) {
        if constexpr (S == sim::physics::bounds::side::left) {
            // move left to right
            double delta;
            for (auto i: indices) {
                auto x = particles[i].getX();
                delta = std::abs(x_0_min - x[0]);
                delta = std::min(r_cutoff, delta);
                x[0] = x_0_max - delta;
                particles[i].setX(x);
                cells[xToCellCoords(i)].emplace_back(i);

            }
        } else if constexpr (S == sim::physics::bounds::side::right) {
            // move right to left
            double delta;
            for (auto i: indices) {
                auto x = particles[i].getX();
                delta = std::abs(x[0] - x_0_max);
                delta = std::min(r_cutoff, delta);
                x[0] = x_0_min + delta;
                particles[i].setX(x);
                cells[xToCellCoords(i)].emplace_back(i);
            }
        } else if constexpr (S == sim::physics::bounds::side::bottom) {
            // move bot to top
            double delta;
            for (auto i: indices) {
                auto x = particles[i].getX();
                delta = std::abs(x_1_min - x[1]);
                delta = std::min(r_cutoff, delta);
                x[1] = x_1_max - delta;
                particles[i].setX(x);
                cells[xToCellCoords(i)].emplace_back(i);
            }
        } else if constexpr (S == sim::physics::bounds::side::top) {
            // move top to bot
            double delta;
            for (auto i: indices) {
                auto x = particles[i].getX();
                delta = std::abs(x[1] - x_1_max);
                delta = std::min(r_cutoff, delta);
                x[1] = x_1_min + delta;
                particles[i].setX(x);
                cells[xToCellCoords(i)].emplace_back(i);
            }
        } else if constexpr (S == sim::physics::bounds::side::front) {
            // move front to rear
            double delta;
            for (auto i: indices) {
                auto x = particles[i].getX();
                delta = std::abs(x_2_min - x[2]);
                delta = std::min(r_cutoff, delta);
                x[2] = x_2_max - delta;
                particles[i].setX(x);
                cells[xToCellCoords(i)].emplace_back(i);
            }
        } else if constexpr (S == sim::physics::bounds::side::rear) {
            // move rear to front
            double delta;
            for (auto i: indices) {
                auto x = particles[i].getX();
                delta = std::abs(x[2] - x_2_max);
                delta = std::min(r_cutoff, delta);
                x[2] = x_2_min + delta;
                particles[i].setX(x);
                cells[xToCellCoords(i)].emplace_back(i);
            }
        }
    }

private:
    /**
     * Computes cell coordinates from particle coordinates
     * index i of x coords
     * */
    unsigned int xToCellCoords(unsigned int i) {
        std::array<unsigned int, 3> cellCoordinate = {0, 0, 0};
        const auto x = particles[i].getX();
        if (x[0] > 0) cellCoordinate[0] = (unsigned int) (x[0] / r_cutoff);
        if (x[1] > 0) cellCoordinate[1] = (unsigned int) (x[1] / r_cutoff);
        if (x[2] > 0) cellCoordinate[2] = (unsigned int) (x[2] / r_cutoff);
        return cellIndexFromCellCoordinates(cellCoordinate);
    }

public:

    /**
     * Stores all particles, that are within the threshold of sixth root of 2 times sigma, into the halo at the opposing side.
     * */
    template<sim::physics::bounds::side S>
    void storeBorderParticlesToHalo(bool mMinor, bool mMajor) {
        if constexpr (S == sim::physics::bounds::side::left) {
            // left x = x_0 = min
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_left = cells[cellIndexFromCellCoordinates({0, x_1, x_2})];
                    for (auto i: cell_indices_left) {
                        unsigned long g0 = gridDimensions[0] + 1;
                        unsigned long g1 = x_1 + 1
                                - mMinor * (x_1 == gridDimensions[1] - 1) * static_cast<long>(gridDimensions[1])
                                - mMinor * (x_1 == gridDimensions[1] - 2) * (static_cast<long>(gridDimensions[1]) - 1)
                                + mMinor * (x_1 == 0) * static_cast<long>(gridDimensions[1]);
                        unsigned long g2 = x_2 + 1
                                - mMajor * (x_2 == gridDimensions[2] - 1) * static_cast<long>(gridDimensions[2])
                                - mMajor * (x_2 == gridDimensions[2] - 2) * (static_cast<long>(gridDimensions[2]) - 1)
                                + mMajor * (x_2 == 0) * static_cast<long>(gridDimensions[2]);
                        cells.store(g0, g1, g2, i);
                        cells.store(g0, x_1 + 1, x_2 + 1, i);
                    }
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::right) {
            // right x = x_0 = max and left x = x_0 = min
            for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_right = cells[cellIndexFromCellCoordinates({gridDimensions[0] - 1, x_1, x_2})];
                    for (auto i: cell_indices_right) {
                        unsigned long g0 = 0;
                        unsigned long g1 = x_1 + 1
                                           - mMinor * (x_1 == gridDimensions[1] - 1) * static_cast<long>(gridDimensions[1])
                                           - mMinor * (x_1 == gridDimensions[1] - 2) * (static_cast<long>(gridDimensions[1]) - 1)
                                           + mMinor * (x_1 == 0) * static_cast<long>(gridDimensions[1]);
                        unsigned long g2 = x_2 + 1
                                           - mMajor * (x_2 == gridDimensions[2] - 1) * static_cast<long>(gridDimensions[2])
                                           - mMajor * (x_2 == gridDimensions[2] - 1) * (static_cast<long>(gridDimensions[2]) - 1)
                                           + mMajor * (x_2 == 0) * static_cast<long>(gridDimensions[2]);
                        cells.store(g0, g1, g2, i);
                        cells.store(g0, x_1 + 1, x_2 + 1, i);
                    }
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::bottom) {
            // top y = x_1 = max and bottom y = x_1 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_bot = cells[cellIndexFromCellCoordinates({x_0, 0, x_2})];
                    for (auto i: cell_indices_bot) {
                        unsigned long g0 = x_0 + 1
                                           - mMinor * (x_0 == gridDimensions[0] - 1) * static_cast<long>(gridDimensions[0])
                                           - mMinor * (x_0 == gridDimensions[0] - 2) * (static_cast<long>(gridDimensions[0]) - 1)
                                           + mMinor * (x_0 == 0) * static_cast<long>(gridDimensions[0]);
                        unsigned long g1 = gridDimensions[1] + 1;
                        unsigned long g2 = x_2 + 1
                                           - mMajor * (x_2 == gridDimensions[2] - 1) * static_cast<long>(gridDimensions[2])
                                           - mMajor * (x_2 == gridDimensions[2] - 2) * (static_cast<long>(gridDimensions[2]) - 1)
                                           + mMajor * (x_2 == 0) * static_cast<long>(gridDimensions[2]);
                        cells.store(g0, g1, g2, i);
                        cells.store(x_0 + 1, g1, x_2 + 1, i);
                    }
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::top) {
            // top y = x_1 = max and bottom y = x_1 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_2 = 0; x_2 < gridDimensions[2]; x_2++) {
                    auto &cell_indices_top = cells[cellIndexFromCellCoordinates({x_0, gridDimensions[1] - 1, x_2})];
                    for (auto i: cell_indices_top) {
                        unsigned long g0 = x_0 + 1
                                           - mMinor * (x_0 == gridDimensions[0] - 1) * static_cast<long>(gridDimensions[0])
                                           - mMinor * (x_0 == gridDimensions[0] - 2) * (static_cast<long>(gridDimensions[0]) - 1)
                                           + mMinor * (x_0 == 0) * static_cast<long>(gridDimensions[0]);
                        unsigned long g1 = 0;
                        unsigned long g2 = x_2 + 1
                                           - mMajor * (x_2 == gridDimensions[2] - 1) * static_cast<long>(gridDimensions[2])
                                           - mMajor * (x_2 == gridDimensions[2] - 2) * (static_cast<long>(gridDimensions[2]) - 1)
                                           + mMajor * (x_2 == 0) * static_cast<long>(gridDimensions[2]);
                        cells.store(g0, g1, g2, i);
                        cells.store(x_0 + 1, g1, x_2 + 1, i);
                    }
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::front) {
            // back z = x_2 = max and front z = x_2 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                    auto &cell_indices_front = cells[cellIndexFromCellCoordinates({x_0, x_1, 0})];
                    for (auto i: cell_indices_front) {
                        unsigned long g0 = x_0 + 1
                                           - mMinor * (x_0 == gridDimensions[0] - 1) * static_cast<long>(gridDimensions[0])
                                           - mMinor * (x_0 == gridDimensions[0] - 2) * (static_cast<long>(gridDimensions[0]) - 1)
                                           + mMinor * (x_0 == 0) * static_cast<long>(gridDimensions[0]);
                        unsigned long g1 = x_1 + 1
                                           - mMajor * (x_1 == gridDimensions[1] - 1) * static_cast<long>(gridDimensions[1])
                                           - mMajor * (x_1 == gridDimensions[1] - 2) * (static_cast<long>(gridDimensions[1]) - 1)
                                           + mMajor * (x_1 == 0) * static_cast<long>(gridDimensions[1]);
                        unsigned long g2 = gridDimensions[2] + 1;
                        cells.store(g0, g1, g2, i);
                        cells.store(x_0 + 1, x_1 + 1, g2, i);
                    }
                }
            }
        } else if constexpr (S == sim::physics::bounds::side::rear) {
            // back z = x_2 = max and front z = x_2 = min
            for (unsigned int x_0 = 0; x_0 < gridDimensions[0]; x_0++) {
                for (unsigned int x_1 = 0; x_1 < gridDimensions[1]; x_1++) {
                    auto &cell_indices_back = cells[cellIndexFromCellCoordinates({x_0, x_1, gridDimensions[2] - 1})];
                    for (auto i: cell_indices_back) {
                        unsigned long g0 = x_0 + 1
                                           - mMinor * (x_0 == gridDimensions[0] - 1) * static_cast<long>(gridDimensions[0])
                                           - mMinor * (x_0 == gridDimensions[0] - 2) * (static_cast<long>(gridDimensions[0]) - 1)
                                           + mMinor * (x_0 == 0) * static_cast<long>(gridDimensions[0]);
                        unsigned long g1 = x_1 + 1
                                           - mMajor * (x_1 == gridDimensions[1] - 1) * static_cast<long>(gridDimensions[1])
                                           - mMajor * (x_1 == gridDimensions[1] - 2) * (static_cast<long>(gridDimensions[1]) - 1)
                                           + mMajor * (x_1 == 0) * static_cast<long>(gridDimensions[1]);
                        unsigned long g2 = 0;
                        cells.store(g0, g1, g2, i);
                        cells.store(x_0 + 1, x_1 + 1, g2, i);
                    }
                }
            }
        }
    }

    /**
     * Clears all halo entries
     * */
    void clearHalo() {
        cells.clearOuter();
    }

    /**
     * @brief getter for gridDimensions.
     * There are gridDimensions[0]*gridDimensions[1]*gridDimensions[2] cells used
     *
     */
    std::array<unsigned int, 3> getGridDimensions();

    /**
     * Runs the function on the internal data
     * */
    template<typename F>
    void runOnData(F fun) {
        fun(particles, membranes, cells, count, activeParticles, id_to_index);
    }

    /**
     * @brief Applies the given function to all Particles
     *
     * @param function
     */
     template<typename F>
    void forAllParticles(F fun) {
        for (unsigned long id: activeParticles) {
            fun(particles[id_to_index[id]]);
        }
    }

    /**
     * @brief Applies given function to all pairs of Particles p_i, p_j, where p_i < p_j once
     *  (If f(p_i, p_j) got invoked, f(p_j, p_i) won't get invoked with the same i and j)
     * @param function
     */
     template <typename F>
    void forAllPairs(F fun) {
        for (u_int32_t i = 0; i < activeSize(); i++) {
            for (u_int32_t j = i + 1; j < activeSize(); j++) {
                fun(particles[id_to_index[activeParticles[i]]], particles[id_to_index[activeParticles[j]]]);
            }
        }
    }

    /**
     * @brief Applies given function to all pairs of Particles that are connected by a spring due to Membranes
     * Also iterates over all Membranes created
     * analogon to forAllPairs
     * This method is too slow for actual use but might be good for debugging purposes or prototyping
     */
    [[maybe_unused]] void forAllMembraneSprings(const std::function<void(Particle &p1, Particle &p2, double desiredDistance, double springStrength)> &function);

    /**
     * Handles interactions between halo and border cells.
     * Runs fun on ever particle pair with one particle being in the halo and the other in the border.
     * Clears handled halo cell
     * */
    template<sim::physics::bounds::side S, typename F>
    void forAllPairsHaloSide(F fun) {
        using p = std::array<unsigned int, 3>;
        using o = std::array<int, 3>;
        //     6 -- 7
        //  2 -- 3  |
        //  |  4 |  5
        //  0 -- 1
        std::array<p, 8> hcCoords = {p{0, 0, 0}, p{gridDimensions[0] + 1, 0, 0}, p{0, gridDimensions[1] + 1, 0}, p{gridDimensions[0] + 1, gridDimensions[1] + 1, 0},
                                     p{0,0,gridDimensions[2]+1}, p{gridDimensions[0]+1, 0, gridDimensions[2]+1}, p{0, gridDimensions[1]+1, gridDimensions[2]+1}, p{gridDimensions[0]+1, gridDimensions[1]+1, gridDimensions[2]+1}};
        constexpr std::array<o, 8> bcOffset = {o{1,1,1}, o{-1,1,1}, o{1,-1,1},o{-1,-1,1}, o{1,1,-1}, o{-1,1,-1}, o{1,-1,-1},o{-1,-1,-1}};
        // index into sideCIndices is S
        using i3 = std::array<int, 3>;
        using i4 = std::array<int, 4>;
        constexpr std::array<i4, 6> sideCIndices = {i4{0, 2, 4, 6}, i4{1,3,5,7}, i4{2,3,6,7}, i4{0,1,4,5}, i4{0,1,2,3}, i4{4,5,6,7}};

        // corner
        {
            for(auto cIndex : sideCIndices[S]) {
                auto& hCell = cells.getOuter(hcCoords[cIndex][0], hcCoords[cIndex][1], hcCoords[cIndex][2]);
                for (auto indexI : hCell) {
                    //transform h coords
                    const auto x = particles[indexI].getX();
                    particles[indexI].add_to_X({
                        (std::min(hcCoords[cIndex][0], 1u)) * domainSize[0] - (1 - std::min(hcCoords[cIndex][0], 1u)) * domainSize[0],
                        (std::min(hcCoords[cIndex][1], 1u)) * domainSize[1] - (1 - std::min(hcCoords[cIndex][1], 1u)) * domainSize[1],
                        (std::min(hcCoords[cIndex][2], 1u)) * domainSize[2] - (1 - std::min(hcCoords[cIndex][2], 1u)) * domainSize[2]
                    });

                    for(int scale0{1}; scale0 < 3; scale0++) {
                        for(int scale1{1}; scale1 < 3; scale1++) {
                            for(int scale2{1}; scale2 < 3; scale2++) {
                                auto& bCell = cells.getOuter(hcCoords[cIndex][0]+bcOffset[cIndex][0]*scale0, hcCoords[cIndex][1]+bcOffset[cIndex][1]*scale1, hcCoords[cIndex][2]+bcOffset[cIndex][2]*scale2);
                                //apply function
                                for (auto indexJ : bCell) {
                                    fun(particles[indexI], particles[indexJ]);
                                }
                            }
                        }
                    }

                    //write back original value
                    particles[indexI].setX(x);
                }
                hCell.clear();
            }
        }
        // edge
        {
            using d = std::array<i3, 3>;
            // tuple<starting_point, dimension, check_directions>
            using t = const std::tuple<std::array<unsigned int, 3>, std::array<unsigned int, 3>, d>;
            const std::array<t, 12> edges {
                    //left to right
                    t{{1,0,0},{gridDimensions[0],1,1},d{i3{0,1,1}, i3{-1,1,1}, i3{1,1,1}}},
                    t{{1,gridDimensions[1]+1,0},{gridDimensions[0],1,1},d{i3{0,-1,1},{-1,-1,1},{1,-1,1}}},
                    t{{1,0,gridDimensions[2]+1},{gridDimensions[0],1,1},d{i3{0,1,-1},{-1,1,-1},{1,1,-1}}},
                    t{{1,gridDimensions[1]+1,gridDimensions[2]+1},{gridDimensions[0],1,1},d{i3{0,-1,-1},{-1,-1,-1},{1,-1,-1}}},
                    //bottom to top
                    t{{0,1,0},{1,gridDimensions[1],1},d{i3{1,0,1},{1,-1,1},{1,1,1}}},
                    t{{gridDimensions[0]+1,1,0},{1,gridDimensions[1],1},d{i3{-1,0,1},{-1,-1,1},{-1,1,1}}},
                    t{{0,1,gridDimensions[2]+1},{1,gridDimensions[1],1},d{i3{1,0,-1},{1,-1,-1},{1,1,-1}}},
                    t{{gridDimensions[0]+1,1,gridDimensions[2]+1},{1,gridDimensions[1],1},d{i3{-1,0,-1},{-1,-1,-1},{-1,1,-1}}},
                    //front to rear
                    t{{0,0,1},{1,1,gridDimensions[2]},d{i3{1,1,0},{1,1,-1},{1,1,1}}},
                    t{{gridDimensions[0]+1,0,1},{1,1,gridDimensions[2]},d{i3{-1,1,0},{-1,1,-1},{-1,1,1}}},
                    t{{0,gridDimensions[1]+1,1},{1,1,gridDimensions[2]},d{i3{1,-1,0},{1,-1,-1},{1,-1,1}}},
                    t{{gridDimensions[0]+1,gridDimensions[1]+1,1},{1,1,gridDimensions[2]},d{i3{-1,-1,0},{-1,-1,-1},{-1,-1,1}}}
            };
            constexpr std::array<i4, 6> sideEIndices = {i4{4,6,8,10}, {5,7,9,11}, {1,3,10,11}, {0,2,8,9}, {0,1,4,5}, {2,3,6,7}};
            for(auto indexE : sideEIndices[S]) {
                auto& [start, dim, dirs] = edges[indexE];

                for(unsigned long x_0 {start[0]}; x_0 < start[0] + dim[0]; x_0++) {
                    for (unsigned long x_1{start[1]}; x_1 < start[1] + dim[1]; x_1++) {
                        for (unsigned long x_2{start[2]}; x_2 < start[2] + dim[2]; x_2++) {
                            auto &hCell = cells.getOuter(x_0, x_1, x_2);

                            for (auto indexI : hCell) {
                                //transform h coords
                                const auto x = particles[indexI].getX();
                                particles[indexI].add_to_X({
                                                                   -domainSize[0] * dirs[0][0],
                                                                   -domainSize[1] * dirs[0][1],
                                                                   -domainSize[2] * dirs[0][2]
                                                           });

                                //go in check direction
                                for (auto &dir: dirs) {
                                    for(int scale0{0}; scale0 < 2; scale0++) {
                                        for(int scale1{0}; scale1 < 2; scale1++) {
                                            for(int scale2{0}; scale2 < 2; scale2++) {
                                                auto &bCell = cells.getInnerGlobal(x_0 + dir[0] + dirs[0][0]*scale0, x_1 + dir[1] + dirs[0][1]*scale1, x_2 + dir[2] + dirs[0][2]*scale2);
                                                for (auto indexJ: bCell) {
                                                    fun(particles[indexI], particles[indexJ]);
                                                }
                                            }
                                        }
                                    }
                                }

                                //write back original value
                                particles[indexI].setX(x);
                            }

                            hCell.clear();
                        }
                    }
                }
            }
        }
        //handle planes
        {
            // tuple<starting_point, dimension, check_direction>
            using ma = std::array<i3, 9>;
            using t = const std::tuple<std::array<unsigned int, 3>, std::array<unsigned int, 3>, ma>;
            const std::array<t, 6> planes {
                    //left
                    t{{0,1,1},{1,gridDimensions[1],gridDimensions[2]},ma{i3{1,1,1}, {1,1,0}, {1,1,-1},
                                                                        i3{1,0,1}, {1,0,0}, {1,0,-1},
                                                                        i3{1,-1,1},{1,-1,0},{1,-1,-1}}},
                    //right
                    t{{gridDimensions[0]+1,1,1},{1,gridDimensions[1],gridDimensions[2]},ma{i3{-1,1,1}, {-1,1,0}, {-1,1,-1},
                                                                                          i3{-1,0,1}, {-1,0,0}, {-1,0,-1},
                                                                                          i3{-1,-1,1},{-1,-1,0},{-1,-1,-1}}},
                    //top
                    t{{1,gridDimensions[1]+1,1},{gridDimensions[0],1,gridDimensions[2]},ma{i3{-1,-1,1}, {0,-1,1}, {1,-1,1},
                                                                                          i3{-1,-1,0}, {0,-1,0}, {1,-1,0},
                                                                                          i3{-1,-1,-1},{0,-1,-1},{1,-1,-1}}},
                    //bottom
                    t{{1,0,1},{gridDimensions[0],1,gridDimensions[2]},ma{i3{-1,1,1}, {0,1,1}, {1,1,1},
                                                                        i3{-1,1,0}, {0,1,0}, {1,1,0},
                                                                        i3{-1,1,-1},{0,1,-1},{1,1,-1}}},
                    //front
                    t{{1,1,0},{gridDimensions[0],gridDimensions[1],1},ma{i3{-1,1,1}, {0,1,1}, {1,1,1},
                                                                        i3{-1,0,1}, {0,0,1}, {1,0,1},
                                                                        i3{-1,-1,1},{0,-1,1},{1,-1,1}}},
                    //rear
                    t{{1,1,gridDimensions[2]+1},{gridDimensions[0],gridDimensions[1],1},ma{i3{-1,1,-1}, {0,1,-1}, {1,1,-1},
                                                                                          i3{-1,0,-1}, {0,0,-1}, {1,0,-1},
                                                                                          i3{-1,-1,-1},{0,-1,-1},{1,-1,-1}}}
            };

            auto& [start, dim, dirs] = planes[S];
            for(unsigned long x_0 {start[0]}; x_0 < start[0] + dim[0]; x_0++) {
                for(unsigned long x_1 {start[1]}; x_1 < start[1] + dim[1]; x_1++){
                    for(unsigned long x_2 {start[2]}; x_2 < start[2] + dim[2]; x_2++){
                        auto& hCell = cells.getOuter(x_0, x_1, x_2);
                            for (auto indexI : hCell) {
                                //transform h coords
                                const auto x = particles[indexI].getX();
                                particles[indexI].add_to_X({
                                                                   -domainSize[0] * dirs[4][0],
                                                                   -domainSize[1] * dirs[4][1],
                                                                   -domainSize[2] * dirs[4][2]
                                                           });

                                //go in check direction
                                for (int scaleOffset{0}; scaleOffset < 2; scaleOffset++) {
                                    for (auto &dir: dirs) {
                                        auto &bCell = cells.getInnerGlobal(x_0 + dir[0] + scaleOffset * dirs[4][0], x_1 + dir[1] + scaleOffset * dirs[4][1], x_2 + dir[2] + scaleOffset * dirs[4][2]);
                                        for (auto indexJ: bCell) {
                                            fun(particles[indexI], particles[indexJ]);
                                        }
                                    }
                                }

                                //write back original value
                                particles[indexI].setX(x);
                            }

                            hCell.clear();
                    }
                }
            }
        }//handle planes
    }

    /**
     * Performs fun on provided data. All lambda args particle container internal data.
     * Will be applied on every cell.
     * */
    template<typename F>
    void forAllCells(F fun) {
        for (auto &cellItems: cells) {
            fun(particles, count, cellItems);
        }
    }

    /**
     * @brief Performs function on all Pairs that are in the same cell
     * 
     * @param function 
     */
    void forAllPairsInSameCell(const std::function<void(Particle &p1, Particle &p2)> &function);

    /**
     * Initializes generateDistinctCellNeighbours cache.
     * */
    void initTaskModel();

    void initAlternativeTaskModel();
private:
    std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>> taskModelCache;
    std::vector<std::pair<unsigned long, unsigned long>> alternativeTaskModelCache;

public:

    /**
     * Generates all cell pairs of neighbours.
     * Pairs that do not interfere with each other are in their own vector.
     * Pairs contain indices to address the cells.
     * vector of task groups -> vector of
     * each task group needs a barrier -> vector of
     * task group is vector task for one thread -> vector of
     * task is vector of pairs -> pair
     * */
    const std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>>& generateDistinctCellNeighbours();

    const std::vector<std::pair<unsigned long, unsigned long>>& generateDistinctAlternativeCellNeighbours();

    /**
     * Performs fun on provided data. All lambda args particle container internal data.
     * Will be applied on every distinct cell neighbours. (Set-Wise) I.e. {a,b} = {b,a}.
     * Arguments:
     *
     * std::vector<double> &force,
     * std::vector<double> &oldForce,
     * std::vector<double> &x,
     * std::vector<double> &v,
     * std::vector<double> &m,
     * std::vector<int> &type,
     * unsigned long count,
     * std::vector<unsigned long> &cell0Items
     * std::vector<unsigned long> &cell1Items
     * std::vector<double> &eps,
     * std::vector<double> &sig
     * */
    template<typename F>
    void forAllDistinctCellNeighbours(F fun) {
        //Implementation3:

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

            for(unsigned int x0 = lowerBounds[c][0]; x0 < upperBounds[c][0]; x0++){
                for(unsigned int x1 = lowerBounds[c][1]; x1 < upperBounds[c][1]; x1++){
                    for(unsigned int x2 = lowerBounds[c][2]; x2 < upperBounds[c][2]; x2++){
                        fun(particles, count,
                            cells[cellIndexFromCellCoordinatesFast(x0, x1, x2)],
                            cells[cellIndexFromCellCoordinatesFast(x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2])]);
                        SPDLOG_TRACE("Cell ({} {} {}) interacted with ({} {} {})", x0, x1, x2, x0 + offsets[c][0], x1 + offsets[c][1], x2 + offsets[c][2]);
                    }
                }
            }

        }
    }

};