#include "sim/analysis/Profiler.h"
#include "data/ParticleContainer.h"
#include "data/Particle.h"

#include <gtest/gtest.h>
#include <vector>
#include <Eigen>

TEST(Profiler, testProfile){
    std::vector<Particle> buffer{};
    buffer.emplace_back(Eigen::Vector3d{5.5,0,0}, Eigen::Vector3d{0,0,0}, 1, 0, 0);
    buffer.emplace_back(Eigen::Vector3d{7.5,0,0}, Eigen::Vector3d{0,-2,0}, 1, 0, 1);
    ParticleContainer pc(buffer, std::array<double, 3>{10., 10., 10.}, 3., {});
    pc.updateCells();
    sim::analysis::Profiler prof{};

    auto profData = prof.profile(10, {10.,10.,10.}, pc);
    ASSERT_EQ(profData.bins, 10);

    for(int i=0; i<10; i++){
        if(i==5){
            ASSERT_EQ(profData.velocities[3*i], 0.);
            ASSERT_EQ(profData.velocities[3*i+1], 0.);
            ASSERT_EQ(profData.velocities[3*i+2], 0.);
        }else if(i==7){
            ASSERT_EQ(profData.velocities[3*i], 0.);
            ASSERT_EQ(profData.velocities[3*i+1], -2.);
            ASSERT_EQ(profData.velocities[3*i+2], 0.);
        }else{
            ASSERT_EQ(profData.densities[i], 0);
            ASSERT_TRUE(std::isnan(profData.velocities[i]));
        }
    }
}