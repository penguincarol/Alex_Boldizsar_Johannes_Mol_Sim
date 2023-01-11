#include "data/ParticleGenerator.h"
#include "data/Membrane.h"
#include "data/Particle.h"
#include "data/Body.h"

#include <gtest/gtest.h>
#include <vector>
#include <Eigen>

/**
 * Test if the Membrane gets initialized as expected even if the buffer is not empty
 */
TEST(ParticleGenerator, generateMembrane){
    std::cout<<"I am here\n";
    std::list<Particle> buf{};
    std::list<Membrane> membrBuf{};


    //fill buffer with some dummy values already
    Eigen::Vector3d x{0,0,0};
    ParticleGenerator::generateParticle(x, x, 2, buf, 0.1, 0.2);
    x[0]+=1;
    ParticleGenerator::generateParticle(x, x, 2, buf, 0.1, 0.2);
    x[0]+=1;
    ParticleGenerator::generateParticle(x, x, 2, buf, 0.1, 0.2);
    ASSERT_EQ(buf.size(), 3);

    //actual generateMembrane stuff:
    Body membr1{membrane, Eigen::Vector3d{1.,1.,1.}, Eigen::Vector3d{1,3,3},
                1.2, 2, Eigen::Vector3d{0,0,0}, 1.2, 5};
    ParticleGenerator::generateMembrane(membr1, 0.1, buf, membrBuf, 3, 0.1, 0.2);

    ASSERT_EQ(buf.size(), 3 + 3*3);
    ASSERT_EQ(membrBuf.size(), 1);

    membrBuf.front().printMembrNodeStructure();

}