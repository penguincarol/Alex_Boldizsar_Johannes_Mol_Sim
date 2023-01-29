#include <gtest/gtest.h>
#include <Eigen>
#include <vector>

#include "data/Particle.h"
#include "data/ParticleContainer.h"
#include "data/ParticleGenerator.h"
#include "sim/physics/force/FMembrane.h"
#include "sim/physics/force/FLennardJonesCells.h"
#include "sim/physics/force/FLennardJonesCellsOMP.h"
#include "sim/physics/force/FMembranePull.h"
#include "sim/physics/position/XStoermerVelvetOMP.h"
#include "sim/physics/position/XStoermerVelvet.h"
#include "sim/physics/velocity/VStoermerVelvetOMP.h"
#include "sim/physics/velocity/VStoermerVelvet.h"

static void vectorEqual(const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){
    ASSERT_DOUBLE_EQ(lhs[0], rhs[0])<< "ForceLeft was "  << lhs << "\nbut was expected to be " << rhs << "\n";
    ASSERT_DOUBLE_EQ(lhs[1], rhs[1])<< "ForceLeft was "  << lhs << "\nbut was expected to be " << rhs << "\n";;
    ASSERT_DOUBLE_EQ(lhs[2], rhs[2])<< "ForceLeft was "  << lhs << "\nbut was expected to be " << rhs << "\n";;
    //return true;
    /*for(int i=0; i<3; i++){
        if(lhs[i]-rhs[i] > std::max(lhs[i]*0.000000001, 0.000'000'000'000'001)){
            return false;
        }
    }
    return true;*/
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
    //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[0][1]).getF(), forceLeft)) << "ForceLeft was "  << pc.getParticle(membrNodes[0][1]).getF() << "\nbut was expected to be " << forceLeft << "\n";    //particle to the left
    vectorEqual(pc.getParticle(membrNodes[0][1]).getF(), forceLeft);
    Eigen::Vector3d forceRight{k * (1.-0.5), 0., 0.};
    //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[2][1]).getF(), forceRight)) << "ForceRight was "  << pc.getParticle(membrNodes[2][1]).getF() << "\nbut was expected to be " << forceRight << "\n";    //particle to the right
    vectorEqual(pc.getParticle(membrNodes[2][1]).getF(), forceRight);

    double normAbove = pc.getParticle(membrNodes[1][2]).getF().norm();  //particle Above
    double aboveDistance = std::sqrt(0.5 * 0.5 + 1 * 1);
    Eigen::Vector3d forceAbove{normAbove * 0.5 / aboveDistance, -normAbove * 1 / aboveDistance, 0.};
    ASSERT_DOUBLE_EQ(normAbove, std::abs(k * (aboveDistance - 1)));
    //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[1][2]).getF(), forceAbove)) << "Force on Particle above was "  << pc.getParticle(membrNodes[1][2]).getF() << "\nbut was expected to be " << forceAbove << "\n";
    vectorEqual(pc.getParticle(membrNodes[1][2]).getF(), forceAbove);

    double normLeftDown = pc.getParticle(membrNodes[0][0]).getF().norm();   //particle bottom left
    double diagDistanceLeftDown = std::sqrt(1.5 * 1.5 + 1 * 1);
    Eigen::Vector3d forceLeftDown{normLeftDown * 1.5 / diagDistanceLeftDown,
                              normLeftDown * 1 / diagDistanceLeftDown, 0.};
    ASSERT_DOUBLE_EQ(normLeftDown, k * (diagDistanceLeftDown - std::sqrt(2) * 1));
    //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[0][0]).getF(), forceLeftDown)) << "Force on the bottom left was "  << pc.getParticle(membrNodes[0][0]).getF() << "\nbut was expected to be " << forceLeftDown << "\n";
    vectorEqual(pc.getParticle(membrNodes[0][0]).getF(), forceLeftDown);

    double normRightUp = pc.getParticle(membrNodes[2][2]).getF().norm();
    double diagDistanceRightUp = std::sqrt(0.5 * 0.5 + 1 * 1);
    Eigen::Vector3d forceRightUp{normRightUp * 0.5 / diagDistanceRightUp,
                                 normRightUp * 1 / diagDistanceRightUp, 0.};
    ASSERT_DOUBLE_EQ(normRightUp, std::abs(k * (diagDistanceRightUp - std::sqrt(2) * 1)));
    //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[2][2]).getF(), forceRightUp));
    vectorEqual(pc.getParticle(membrNodes[2][2]).getF(), forceRightUp);

    Eigen::Vector3d fMiddleParticle{-forceLeft[0]-forceRight[0]- 2*forceLeftDown[0] - 2*forceRightUp[0] - 2*forceAbove[0], 0., 0.};
    //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[1][1]).getF(), fMiddleParticle));
    vectorEqual(pc.getParticle(membrNodes[1][1]).getF(), fMiddleParticle);
}


/***
 * This test is supposed to check whether LennardJonesCells is getting truncated for Membranes correctly
 */
TEST(FMembrane, LennardJonesTruncateIfWanted) {
    std::list<Particle> buf;
    std::list<Membrane> membrBuf;

    //desiredDistance = startingDistance
    Body membr;
    membr.shape = membrane;
    membr.fixpoint = {6, 1, 1};
    membr.dimensions = {4,4, 1};
    membr.distance = 3.05;
    membr.mass = 1;
    membr.start_velocity = {0., 0., 0.};
    membr.desiredDistance = 3.05;
    membr.springStrength = 3;
    membr.pullEndTime = 0;
    membr.pullForce = {0, 0, 0};
    membr.pullIndices = {};

    double rt2_6 = std::pow(2.0, 1.0/6.0);
    double sig{3.0/rt2_6};      //=> truncation should happen at 3 since std::pow(20, 1.0/6.0)*sig = 3
    double eps{sig*1.2};

    ParticleGenerator::generateMembrane(membr, 0, buf, membrBuf, 3, sig, eps);

    std::vector bufVec(buf.begin(), buf.end());
    std::vector<Membrane> membrVec(membrBuf.begin(), membrBuf.end());
    ParticleContainer pc(bufVec, std::array<double, 3>{10., 10., 10.}, 3., membrVec);

    //auto fMem = sim::physics::force::FMembrane(0, 100, 0.01, eps, sig, pc);
    //intentionally giving them absurdly high default eps and sig so that they need to use the membrane eps and sig to truncate
    auto fLen = sim::physics::force::FLennardJonesCells(0, 100, 0.01, 500, 500, pc);
    auto fLenOMP = sim::physics::force::FLennardJonesCells(0, 100, 0.01, 500, 500, pc);

    //since membrane got initialized with 3.05 distance and truncation should take place at a distance of 3
    //no forces should be applied in either case here
    fLen.operator()();
    pc.forAllParticles([&](Particle &p) {
        ASSERT_EQ(p.getF().norm(), 0.);
    });
    fLenOMP.operator()();
    pc.forAllParticles([&](Particle &p) {
        ASSERT_EQ(p.getF().norm(), 0.);
    });
}

/**
 * This test is a sanity check that looks whether you are NOT truncating
 * if you aren't supposed to truncate
 */
TEST(FMembrane, LennardJonesDontAlwaysTruncate) {
    std::list<Particle> buf;
    std::list<Membrane> membrBuf;

    //desiredDistance = startingDistance
    Body membr;
    membr.shape = membrane;
    membr.fixpoint = {6, 1, 1};
    membr.dimensions = {4,4, 1};
    membr.distance = 3.0;
    const double d = membr.distance;
    membr.mass = 1;
    membr.start_velocity = {0., 0., 0.};
    membr.desiredDistance = 3.0;
    membr.springStrength = 5;
    membr.pullEndTime = 0;
    membr.pullForce = {0, 0, 0};
    membr.pullIndices = {};

    double rt2_6 = std::pow(2.0, 1.0/6.0);
    double sigFact=1.00025;
    double sig{(3.0/rt2_6)*sigFact};
    double eps{4.0};

    ParticleGenerator::generateMembrane(membr, 0, buf, membrBuf, 3, sig, eps);

    std::vector bufVec(buf.begin(), buf.end());
    std::vector<Membrane> membrVec(membrBuf.begin(), membrBuf.end());
    ParticleContainer pc(bufVec, std::array<double, 3>{10., 10., 10.}, 3., membrVec);

    //auto fMem = sim::physics::force::FMembrane(0, 100, 0.01, eps, sig, pc);
    //setting eps and sig so that he truncates if he accidentally uses them instead
    auto fLen = sim::physics::force::FLennardJonesCells(0, 100, 0.01, 0.1, 0.1, pc);
    auto fLenOMP = sim::physics::force::FLennardJonesCells(0, 100, 0.01, 0.1, 0.1, pc);

    fLen.operator()();
    //check one outer particle and one inner particle as sample if it works properly
    //auto temp = pc.getMembranes()[0].getMembrNodes()[0][0];
    auto idOut = pc.getMembranes()[0].getMembrNodes()[1][0];
    auto idIn = pc.getMembranes()[0].getMembrNodes()[1][1];
    //Eigen::Vector3d forceOnOut = {0., -(24.0*eps/d*d)* d * ((1.*std::pow(sigFact, 6.0)/2.) - 2*(std::pow(sigFact, 12.0)/4.0)), 0.};
    Eigen::Vector3d forceOnOut = {0., -(24.0*eps/d)* d * ((1.*std::pow(sigFact, 6.0)/2.) - 2*(std::pow(sigFact, 12.0)/4.0)), 0.};

    //ASSERT_TRUE(vectorEqual(pc.getParticle(idIn).getF(), {0,0,0}));
    //ASSERT_TRUE(vectorEqual(pc.getParticle(idOut).getF(), forceOnOut));
    vectorEqual(pc.getParticle(idIn).getF(), {0,0,0});
    vectorEqual(pc.getParticle(idOut).getF(), forceOnOut);

    pc.forAllParticles([&](Particle &p) {
        p.setF({0,0,0});
    });

    fLenOMP.operator()();
    //ASSERT_TRUE(vectorEqual(pc.getParticle(idIn).getF(), {0,0,0}));
    //ASSERT_TRUE(vectorEqual(pc.getParticle(idOut).getF(), forceOnOut));
    vectorEqual(pc.getParticle(idIn).getF(), {0,0,0});
    vectorEqual(pc.getParticle(idOut).getF(), forceOnOut);
}

/**
 * Rounding errors in starting position escalated into rapid oscillations that broke the sim
 * This test is supposed to analyse that behaviour
 */
TEST(FMembrane, increasingErrors){
    std::list<Particle> buf;
    std::list<Membrane> membrBuf;

    //desiredDistance = startingDistance
    Body membr;
    membr.shape = membrane;
    membr.fixpoint = {1, 1, 1};
    membr.dimensions = {2,1, 1};
    membr.distance = 2.2;
    const double d = membr.distance;
    membr.mass = 1;
    membr.start_velocity = {0., 0., 0.};
    membr.desiredDistance = 2.2;
    membr.springStrength = 300;
    membr.pullEndTime = 0;
    membr.pullForce = {0, 0, 0};
    membr.pullIndices = {};

    ParticleGenerator::generateMembrane(membr, 0, buf, membrBuf, 3, 1.2, 1.2);

    std::vector bufVec(buf.begin(), buf.end());
    std::vector<Membrane> membrVec(membrBuf.begin(), membrBuf.end());
    ParticleContainer pc(bufVec, std::array<double, 3>{100., 100., 100.}, 3., membrVec);

    auto fMem = sim::physics::force::FMembrane(0, 100, 0.01, 1., 1., pc);
    auto fLen = sim::physics::force::FLennardJonesCells(0, 100, 0.01, 1., 1., pc);

    //auto xCalc = sim::physics::position::XStoermerVelvetOMP(0., 100, 0.01, 1, 1, pc);
    //auto vCalc = sim::physics::velocity::VStoermerVelvetOMP(0., 100., 0.01, 1, 1, pc);
    auto xCalc = sim::physics::position::XStoermerVelvet(0., 100, 0.01, 1., 1., pc);
    auto vCalc = sim::physics::velocity::VStoermerVelvet(0., 100., 0.01, 1., 1., pc);

    fMem.operator()();
    //ASSERT_DOUBLE_EQ(pc.getParticle(pc.getMembranes()[0].getMembrNodes()[1][0]).getF().norm(), 0.)<<"1 force calc";

    for(size_t i{0}; i < 5'000; i++){
        xCalc.operator()();
        auto dist=(pc.getParticle(pc.getMembranes()[0].getMembrNodes()[1][0]).getX() - pc.getParticle(pc.getMembranes()[0].getMembrNodes()[0][0]).getX()).norm();
        //ASSERT_TRUE(std::abs(dist-membr.desiredDistance)<=0.10000000000000009)<<"Particles at iteration "<<i<<" were " << dist-membr.desiredDistance << " apart from each other";
        pc.updateCells();
        pc.clearStoreForce();
        fMem.operator()();
        fLen.operator()();
        vCalc.operator()();
    }
    ASSERT_DOUBLE_EQ(pc.getParticle(pc.getMembranes()[0].getMembrNodes()[1][0]).getF().norm(), 0.)<< "Force at " << 1 << " " << 0 << "after 1k iterations unequal to 0";
    ASSERT_DOUBLE_EQ(pc.getParticle(pc.getMembranes()[0].getMembrNodes()[0][0]).getF().norm(), 0.)<< "Force at " << 0 << " " << 0 << "after 1k iterations unequal to 0";
    //ASSERT_DOUBLE_EQ(pc.getParticle(pc.getMembranes()[0].getMembrNodes()[2][2]).getF().norm(), 0.)<< "Force at " << 2 << " " << 2 << "after 1k iterations unequal to 0";
    //ASSERT_DOUBLE_EQ(pc.getParticle(pc.getMembranes()[0].getMembrNodes()[1][1]).getF().norm(), 0.)<< "Force at " << 1 << " " << 1 << "after 1k iterations unequal to 0";
    //ASSERT_DOUBLE_EQ(pc.getParticle(pc.getMembranes()[0].getMembrNodes()[1][0]).getF().norm(), 0.)<< "Force at " << 1 << " " << 0<< "after 1k iterations unequal to 0";

    /*
    for(size_t i{0}; i < pc.getMembranes()[0].getMembrNodes().size(); i++){
        for(size_t j{0}; j < pc.getMembranes()[0].getMembrNodes()[0].size(); j++){
            ASSERT_DOUBLE_EQ(pc.getParticle(pc.getMembranes()[0].getMembrNodes()[i][j]).getF().norm(), 0.)<< "Force at " << i << " " << j << "after 1k iterations unequal to 0";
        }
    }*/
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

    //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[1][1]).getF(), Eigen::Vector3d{0, 1, 2}));
    //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[3][1]).getF(), Eigen::Vector3d{0, 1, 2}));
    vectorEqual(pc.getParticle(membrNodes[1][1]).getF(), Eigen::Vector3d{0, 1, 2});
    vectorEqual(pc.getParticle(membrNodes[3][1]).getF(), Eigen::Vector3d{0, 1, 2});

    for(int i=0; i< 4; i++){
        for(int j=0; j<4; j++){
            if((i!=3 && i!=1) || j!=1){
                //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[i][j]).getF(), Eigen::Vector3d{0, 0, 0}));
                vectorEqual(pc.getParticle(membrNodes[i][j]).getF(), Eigen::Vector3d{0, 0, 0});
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
            //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[i][j]).getF(), Eigen::Vector3d{0, 0, 0}));
            vectorEqual(pc.getParticle(membrNodes[i][j]).getF(), Eigen::Vector3d{0, 0, 0});
        }
    }

    fPull.operator()();

    for(int i=0; i< 4; i++){
        for(int j=0; j<4; j++){
            //ASSERT_TRUE(vectorEqual(pc.getParticle(membrNodes[i][j]).getF(), Eigen::Vector3d{0, 0, 0}));
            vectorEqual(pc.getParticle(membrNodes[i][j]).getF(), Eigen::Vector3d{0, 0, 0});
        }
    }
}