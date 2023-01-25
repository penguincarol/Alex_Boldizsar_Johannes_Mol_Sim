#include <gtest/gtest.h>
#include <Eigen>
#include <vector>

#include "data/Particle.h"
#include "data/ParticleContainer.h"
#include "data/ParticleGenerator.h"
#include "sim/physics/force/FMembrane.h"
#include "sim/physics/force/FMembranePull.h"

static bool vectorEqual(const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){
    for(int i=0; i<3; i++){
        if(lhs[i]-rhs[i] > std::max(lhs[i]*0.000000001, 0.000000001)){
            return false;
        }
    }
    return true;
}


TEST(FMembrane, operator) {
    std::list<Particle> buf;
    std::list<Membrane> membrBuf;

    Body cub;
    cub.shape = cuboid;
    cub.fixpoint = {1, 1, 1};
    cub.dimensions = {4, 4, 4};
    cub.distance = 1;
    cub.mass = 0.5;
    cub.start_velocity = {0., 0., 0.};

    //desiredDistance = startingDistance
    Body membr;
    membr.shape = membrane;
    membr.fixpoint = {6, 1, 1};
    membr.dimensions = {4,4, 1};
    membr.distance = 1;
    membr.mass = 1;
    membr.start_velocity = {0., 0., 0.};
    membr.desiredDistance = 1;
    membr.springStrength = 3;
    const auto k = membr.springStrength;
    membr.pullEndTime = 0;
    membr.pullForce = {0, 0, 0};
    membr.pullIndices = {};

    ParticleGenerator::generateCuboid(cub, 0, buf, 3, 1, 1);
    ParticleGenerator::generateMembrane(membr, 0, buf, membrBuf, 3, 1, 1);

    std::vector bufVec(buf.begin(), buf.end());
    std::vector<Membrane> membrVec(membrBuf.begin(), membrBuf.end());
    ParticleContainer pc(bufVec, std::array<double, 3>{10., 10., 10.}, 3., membrVec);

    auto fMem = sim::physics::force::FMembrane(0, 3, 0.1, 1, 1, pc);

    fMem.operator()();

    //no force should be there
    pc.forAllParticles([&](Particle &p) {
        ASSERT_EQ(p.getF().norm(), 0.);
    });

    auto membrNodes = pc.getMembranes()[0].getMembrNodes();

    //move one particle in membrane 0.5 to the right
    pc.runOnData([&membrNodes](std::vector<double>&, std::vector<double>&, std::vector<double>& x, std::vector<double>&,
                               std::vector<double>&, std::vector<int>&, unsigned long,
                               std::vector<double>&, std::vector<double>&) {
        auto index = membrNodes[1][1];
        x[3*index + 0] += 0.5;
    });

    fMem.operator()();

    Eigen::Vector3d forceLeft{k * (1.5 - 1), 0., 0.};
    ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[0][1]).getF(), forceLeft)) << "ForceLeft was "  << pc.getParticle(membrNodes[0][1]).getF() << "\nbut was expected to be " << forceLeft << "\n";    //particle to the left
    Eigen::Vector3d forceRight{k * (1.-0.5), 0., 0.};
    ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[2][1]).getF(), forceRight)) << "ForceRight was "  << pc.getParticle(membrNodes[2][1]).getF() << "\nbut was expected to be " << forceRight << "\n";    //particle to the right

    double normAbove = pc.getParticle(membrNodes[1][2]).getF().norm();  //particle Above
    double aboveDistance = std::sqrt(0.5 * 0.5 + 1 * 1);
    Eigen::Vector3d forceAbove{normAbove * 0.5 / aboveDistance, -normAbove * 1 / aboveDistance, 0.};
    ASSERT_DOUBLE_EQ(normAbove, std::abs(k * (aboveDistance - 1)));
    ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[1][2]).getF(), forceAbove)) << "Force on Particle above was "  << pc.getParticle(membrNodes[1][2]).getF() << "\nbut was expected to be " << forceAbove << "\n";

    double normLeftDown = pc.getParticle(membrNodes[0][0]).getF().norm();   //particle bottom left
    double diagDistanceLeftDown = std::sqrt(1.5 * 1.5 + 1 * 1);
    Eigen::Vector3d forceLeftDown{normLeftDown * 1.5 / diagDistanceLeftDown,
                              normLeftDown * 1 / diagDistanceLeftDown, 0.};
    ASSERT_DOUBLE_EQ(normLeftDown, k * (diagDistanceLeftDown - std::sqrt(2) * 1));
    ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[0][0]).getF(), forceLeftDown)) << "Force on the bottom left was "  << pc.getParticle(membrNodes[0][0]).getF() << "\nbut was expected to be " << forceLeftDown << "\n";

    double normRightUp = pc.getParticle(membrNodes[2][2]).getF().norm();
    double diagDistanceRightUp = std::sqrt(0.5 * 0.5 + 1 * 1);
    Eigen::Vector3d forceRightUp{normRightUp * 0.5 / diagDistanceRightUp,
                                 normRightUp * 1 / diagDistanceRightUp, 0.};
    ASSERT_DOUBLE_EQ(normRightUp, std::abs(k * (diagDistanceRightUp - std::sqrt(2) * 1)));
    ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[2][2]).getF(), forceRightUp));

    Eigen::Vector3d fMiddleParticle{-forceLeft[0]-forceRight[0]- 2*forceLeftDown[0] - 2*forceRightUp[0] - 2*forceAbove[0], 0., 0.};
    ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[1][1]).getF(), fMiddleParticle));

}

TEST(FMembranePull, operator) {
    std::list<Particle> buf;
    std::list<Membrane> membrBuf;

    Body membr;
    membr.shape = membrane;
    membr.fixpoint = {6, 1, 1};
    membr.dimensions = {4,4, 1};
    membr.distance = 1;
    membr.mass = 1;
    membr.start_velocity = {0., 0., 0.};
    membr.desiredDistance = 1;
    membr.springStrength = 1;
    membr.pullEndTime = 0.5;
    membr.pullForce = {0, 1, 2};
    membr.pullIndices = {std::array<size_t, 2>{1,1}, std::array<size_t, 2>{3,1}};

    ParticleGenerator::generateMembrane(membr, 0, buf, membrBuf, 3, 1, 1);

    std::vector bufVec(buf.begin(), buf.end());
    std::vector<Membrane> membrVec(membrBuf.begin(), membrBuf.end());
    ParticleContainer pc(bufVec, std::array<double, 3>{10., 10., 10.}, 3., membrVec);

    auto fPull = sim::physics::force::FMembranePull(0, 3, 1, 1, 1, pc);

    pc.forAllParticles([&](Particle& p){
        ASSERT_EQ(p.getF().norm(), 0.);
    });

    fPull.operator()();

    auto membrNodes = pc.getMembranes()[0].getMembrNodes();

    ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[1][1]).getF(), Eigen::Vector3d{0, 1, 2}));
    ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[3][1]).getF(), Eigen::Vector3d{0, 1, 2}));

    for(int i=0; i< 4; i++){
        for(int j=0; j<4; j++){
            if((i!=3 && i!=1) || j!=1){
                ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[i][j]).getF(), Eigen::Vector3d{0, 0, 0}));
            }
        }
    }

    pc.runOnData([&](std::vector<double>& force,
                 std::vector<double>& oldForce,
                 std::vector<double>& x,
                 std::vector<double>& v,
                 std::vector<double>& m,
                 std::vector<int>& type,
                 unsigned long,
                 std::vector<double>&,
                 std::vector<double>&){
        std::for_each(force.begin(), force.end(), [&](double& e){e=0;});
    });

    for(int i=0; i< 4; i++){
        for(int j=0; j<4; j++){
            ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[i][j]).getF(), Eigen::Vector3d{0, 0, 0}));
        }
    }

    fPull.operator()();

    for(int i=0; i< 4; i++){
        for(int j=0; j<4; j++){
            ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[i][j]).getF(), Eigen::Vector3d{0, 0, 0}));
        }
    }
}