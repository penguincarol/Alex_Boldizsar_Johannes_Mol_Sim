//
// Created by alex on 14.01.23.
//

#include <fstream>
#include "Profiler.h"

namespace sim::analysis {
    void
    Profiler::run(size_t bins, std::array<double, 3> domainSize, ParticleContainer &pc, const std::string &fileName) {
        std::string tmp{fileName};
        if (!tmp.ends_with(".csv")) tmp.append(".csv");

        writeCSV(tmp, profile(bins, domainSize, pc));
    }

    Profiler::ProfileData Profiler::profile(size_t bins, std::array<double, 3> domainSize, ParticleContainer &pc) {
        std::vector<double> densities;
        std::vector<double> velocities;
        densities.resize(bins);
        velocities.resize(3 * bins);

        const double delta_x0 = domainSize[0] / bins;
        double x0 = 0;
        pc.runOnActiveData([&](auto &, auto &,
                               Kokkos::View<double*> &x, Kokkos::View<double*> &v,
                               auto &, auto &, auto &, auto &, auto &,
                               std::vector<unsigned long> &activeParticles) {
            for (size_t i{0}; i < bins; i++) {
                size_t count = 0;
                size_t size = activeParticles.size();
                for(size_t j = 0; j < size; j++) {
                    size_t index = activeParticles[j];
                    if (x[3 * index + 0] < x0 || x[3 * index + 0] >= x0 + delta_x0) continue;
                    count++;
                    densities[i] += 1;
                    velocities[3 * i + 0] += v[3 * index + 0];
                    velocities[3 * i + 1] += v[3 * index + 1];
                    velocities[3 * i + 2] += v[3 * index + 2];
                }

                densities[i] /= delta_x0 * domainSize[1] * domainSize[2];
                velocities[3 * i + 0] /= static_cast<double>(count);
                velocities[3 * i + 1] /= static_cast<double>(count);
                velocities[3 * i + 2] /= static_cast<double>(count);
                x0 += delta_x0;
            }
        });


        return {bins, std::move(densities), std::move(velocities)};
    }

    void Profiler::writeCSV(const std::string &fileName, const Profiler::ProfileData &&data) {
        std::ofstream file;
        file.open(fileName);
        file << "Bin;Density;Velocity" << std::endl;
        for (size_t i{0}; i < data.bins; i++) {
            file << i << ";" << data.densities[i] << ";" <<
                 "[" << data.velocities[3 * i + 0] << ","
                     << data.velocities[3 * i + 1] << ","
                     << data.velocities[3 * i + 2] << "]" << std::endl;
        }
        file.close();
    }
} // analysis