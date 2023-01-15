//
// Created by alex on 14.01.23.
//

#pragma once

#include <array>
#include "data/ParticleContainer.h"

namespace sim::analysis {

    /**
     * Used to create a density and velocity profile along the x0-axis of the domain of the simulation
     * */
    class Profiler {
    private:
        struct ProfileData {
            size_t bins;
            std::vector<double> densities;
            std::vector<double> velocities;
        };
    public:
        /**
         * @brief Performs the analysis on the given data in pc.
         * @param bins Bin count
         * @param domainSize Simulation Domain Size
         * */
        static ProfileData profile(size_t bins, std::array<double,3> domainSize, ParticleContainer& pc);

        /**
         * Writes the given data into csv format
         * */
        static void writeCSV(const std::string& fileName, const ProfileData&& data);

        /**
         * Runs the analysis and writes data into csv format
         * */
        static void run(size_t bins, std::array<double,3> domainSize, ParticleContainer& pc, const std::string& fileName);
    };

} // analysis

