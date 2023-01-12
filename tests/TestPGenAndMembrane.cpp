#include "data/ParticleGenerator.h"
#include "data/Membrane.h"
#include "data/Particle.h"
#include "data/Body.h"

#include <gtest/gtest.h>
#include <vector>
#include <Eigen>

bool vectorEqual(const Eigen::Vector3d &lhs, const Eigen::Vector3d &rhs){
    return (lhs[0]==rhs[0] && lhs[1]==rhs[1] && lhs[2]==rhs[2]);
}

/**
 * Test if the Membrane gets initialized as expected even if the buffer is not empty
 */
TEST(ParticleGenerator, generateMembrane){
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

    //check sporadically if some positions of particles are correct
    Membrane& membr = membrBuf.front();

    std::vector<Particle> bufVec(buf.begin(), buf.end());
    std::vector<Membrane> bufMembrVec(membrBuf.begin(), membrBuf.end());

    ASSERT_TRUE(vectorEqual(bufVec[bufMembrVec[0].getMembrNodes()[0][0]].getX(), {1.,1.,1}));    //we take the particle behind the index that membr[0][0] points to and see if it has the expected value
    ASSERT_TRUE(vectorEqual(bufVec[bufMembrVec[0].getMembrNodes()[2][2]].getX(), {1. , 1.+ 2*1.2, 1+ 2*1.2}));
    ASSERT_TRUE(vectorEqual(bufVec[bufMembrVec[0].getMembrNodes()[1][0]].getX(), {1. , 1.+ 1*1.2, 1+ 0*1.2}));

}