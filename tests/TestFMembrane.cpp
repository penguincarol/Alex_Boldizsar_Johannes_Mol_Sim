#include <gtest/gtest.h>
#include <Eigen>
#include <vector>

#include "data/Particle.h"
#include "io/output/Logging.h"
#include "data/ParticleContainer.h"
#include "data/ParticleGenerator.h"
#include "sim/physics/force/FMembrane.h"



TEST(FLennardJonesCells, operator) {
    std::list<Particle> buf;
    std::list<Membrane> membrBuf;

    Body cub;
    cub.shape = cuboid; cub.fixpoint = {1,1,1}; cub.dimensions = {4,4,4}; cub.distance = 1; cub.mass = 0.5; cub.start_velocity = {0., 0., 0.};

    //first case: desiredDistance = startingDistance
    Body membr;
    membr.shape = membrane; membr.fixpoint = {6,1,1}; membr.dimensions = {4,1,4}; membr.distance = 1; membr.mass = 1; membr.start_velocity = {0.,0., 0.};
    membr.desiredDistance = 1; membr.springStrength = 1; membr.pullEndTime = 0; membr.pullForce = {0,0,0}; membr.pullIndices = {};

    ParticleGenerator::generateCuboid(cub, 0, buf, 3, 1, 1);
    ParticleGenerator::generateMembrane(membr, 0, buf, membrBuf, 3, 1, 1);

    std::vector bufVec(buf.begin(), buf.end());
    std::vector<Membrane> membrVec(membrBuf.begin(), membrBuf.end());
    ParticleContainer pc(bufVec,  std::array<double, 3> {10.,10.,10.}, 3., membrVec);

    auto fMem = sim::physics::force::FMembrane(0, 3, 0.1, 1, 1, pc);

    fMem.operator()();

    //no force should be there
    pc.forAllParticles([&](Particle& p){
        ASSERT_EQ(p.getF().norm(), 0.);
    });

}