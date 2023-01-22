#include <gtest/gtest.h>
#include <cmath>

#include "sim/physics/thermostat/Thermostat.h"
#include "data/Particle.h"
#include "data/ParticleGenerator.h"
#include "data/Body.h"

TEST(Thermostat, helperFunctionsAndCooling) {
    std::vector<Particle> buffer{};
    std::array<std::array <double, 3>, 3> velocities;
    velocities[0] = {5, 0., 0.};
    velocities[1] = {0.,0.,0.};
    velocities[2] = {0.,5/std::sqrt(2),5/std::sqrt(2)}; //v[2] has length 5 (pointing diagonally)

    Particle p1{std::array<double, 3> {1.,1.,1.}, velocities[0], 4., 0};
    Particle p2{std::array<double, 3> {2.,2.,2.}, velocities[1], 1., 0};
    Particle p3{std::array<double, 3> {2.,1.,1.}, velocities[2], 5., 0};
    p1.setID(0);
    p2.setID(1);
    p3.setID(2);
    buffer.emplace_back(p1);
    buffer.emplace_back(p2);
    buffer.emplace_back(p3);

    ParticleContainer pc(buffer, {4.,4.,4.}, 1.0, {}, false);
    //Thermostat(ParticleContainer& particleContainer, double T_t, 
    //unsigned int cT = 100, unsigned int dimensions = 2, double dT = std::numeric_limits<double>::infinity()):
    Thermostat ts(pc, 22, 2, 3, 2, 25, true);

    //temperature ca. 25
    ASSERT_TRUE(std::abs(ts.computeCurrentTemp() - 25) < 25 * 0.000000001)<<"current temp was " << ts.computeCurrentTemp() << " and not approximately 25" << " active Particles: "<< pc.activeSize()<<"\n";

    //vectors unchanged so far
    for(int i{0}; i < 3; i++){
        EXPECT_LE(std::abs(pc.getParticle(i).getV()[0]-velocities[i][0]), 0.0000001 * velocities[i][0]);
        EXPECT_LE(std::abs(pc.getParticle(i).getV()[1]-velocities[i][1]), 0.0000001 * velocities[i][1]);
        EXPECT_LE(std::abs(pc.getParticle(i).getV()[2]-velocities[i][2]), 0.0000001 * velocities[i][2]);
    }

    ts.notify();

    //vectors unchanged so far
    for(int i{0}; i < 3; i++){
        ASSERT_TRUE(pc.getParticle(i).getV()[0] ==  velocities[i][0] && pc.getParticle(i).getV()[1] ==  velocities[i][1] && pc.getParticle(i).getV()[2] ==  velocities[i][2]);
    }

    ts.notify();


    //temperature ca. 23
    ASSERT_TRUE(std::abs(ts.computeCurrentTemp() - 23) < 23 * 0.000000001) << "current temp was " << ts.computeCurrentTemp() << " instead of 23\n";

    ts.notify();
    ts.notify();

    ASSERT_TRUE(std::abs(ts.computeCurrentTemp() - 22) < 22 * 0.000000001) << "current temp was " << ts.computeCurrentTemp() << " instead of 22\n";

}

TEST(Thermostat, Heating) {
    std::vector<Particle> buffer{};
    std::array<std::array <double, 3>, 3> velocities;
    velocities[0] = {5, 0., 0.};
    velocities[1] = {0.,0.,0.};
    velocities[2] = {0.,5/std::sqrt(2),5/std::sqrt(2)}; //v[2] has length 5 (pointing diagonally)

    Particle p1{std::array<double, 3> {1.,1.,1.}, velocities[0], 4., 0};
    Particle p2{std::array<double, 3> {2.,2.,2.}, velocities[1], 1., 0};
    Particle p3{std::array<double, 3> {2.,1.,1.}, velocities[2], 5., 0};
    p1.setID(0);
    p2.setID(1);
    p3.setID(2);
    buffer.emplace_back(p1);
    buffer.emplace_back(p2);
    buffer.emplace_back(p3);

    ParticleContainer pc(buffer, {4.,4.,4.}, 1.0);
    //Thermostat(ParticleContainer& particleContainer, double T_t,
    //unsigned int cT = 100, unsigned int dimensions = 2, double dT = std::numeric_limits<double>::infinity()):
    Thermostat ts(pc, 30, 2, 3, 3, 25, true);

    ts.notify();
    //vectors unchanged so far
    for(int i{0}; i < 3; i++){
        ASSERT_TRUE(pc.getParticle(i).getV()[0] ==  velocities[i][0] && pc.getParticle(i).getV()[1] ==  velocities[i][1] && pc.getParticle(i).getV()[2] ==  velocities[i][2]);
    }

    ts.notify();

    //temperature ca. 28
    ASSERT_TRUE(std::abs(ts.computeCurrentTemp() - 28) < 28 * 0.000000001) << "current temp was " << ts.computeCurrentTemp() << " instead of 28\n";

    ts.notify();
    ts.notify();

    //temperature ca. 30
    ASSERT_TRUE(std::abs(ts.computeCurrentTemp() - 30) < 30 * 0.000000001) << "current temp was " << ts.computeCurrentTemp() << " instead of 30\n";
}

TEST(Thermostat, deltaTInf) {
    std::vector<Particle> buffer{};
    std::array<std::array <double, 3>, 3> velocities;
    velocities[0] = {5, 0., 0.};
    velocities[1] = {0.,0.,0.};
    velocities[2] = {0.,5/std::sqrt(2),5/std::sqrt(2)}; //v[2] has length 5 (pointing diagonally)

    Particle p1{std::array<double, 3> {1.,1.,1.}, velocities[0], 4., 0};
    Particle p2{std::array<double, 3> {2.,2.,2.}, velocities[1], 1., 0};
    Particle p3{std::array<double, 3> {2.,1.,1.}, velocities[2], 5., 0};
    p1.setID(0);
    p2.setID(1);
    p3.setID(2);
    buffer.emplace_back(p1);
    buffer.emplace_back(p2);
    buffer.emplace_back(p3);

    ParticleContainer pc(buffer, {4.,4.,4.}, 1.0);
    //Thermostat(ParticleContainer& particleContainer, double T_t,
    //unsigned int cT = 100, unsigned int dimensions = 2, double dT = std::numeric_limits<double>::infinity()):
    Thermostat ts(pc, 30, 2, 3, std::numeric_limits<double>::infinity(), 25, true);

    ts.notify();
    ts.notify();

    //temperature ca. 30
    ASSERT_TRUE(std::abs(ts.computeCurrentTemp() - 30) < 30 * 0.000000001) << "current temp was " << ts.computeCurrentTemp() << " instead of 30\n";


    //now nothing should happen:
    ts.notify();
    ts.notify();

    //temperature ca. 30
    ASSERT_TRUE(std::abs(ts.computeCurrentTemp() - 30) < 30 * 0.000000001) << "current temp was " << ts.computeCurrentTemp() << " instead of 30\n";
}

TEST(Thermostat, HeatInit){
    std::vector<Particle> buffer{};
    size_t numParWeight1{100000};
    size_t numParWeight5{100000};
    for(size_t i{0}; i<numParWeight1; i++){
        buffer.emplace_back(Eigen::Vector3d{((double )i)/1000, 0, 0}, Eigen::Vector3d{0.,0.,0.}, 1.0, 0);
    }
    for(size_t i{0}; i<numParWeight5; i++){
        buffer.emplace_back(Eigen::Vector3d{((double )i)/1000, 1, 1}, Eigen::Vector3d{0.,0.,0.}, 4.0, 0);
    }

    ParticleContainer pc(buffer, {100., 100., 100.}, 10.0);
    double TInit{30};
    Thermostat ts(pc, 0, 2, 3, 3, 30, true);

    ASSERT_TRUE(std::abs(ts.computeCurrentTemp()-TInit) < TInit*0.0001)<< "Current temp was " << ts.computeCurrentTemp() <<" instead of being approximately " << TInit<< " after initialization with Thermostat"<<std::endl;
}

TEST(Thermostat, pipeFeature){
    Body pipeWall;
    pipeWall.shape = cuboid;
    pipeWall.fixpoint = Eigen::Vector3d{0,0,0};
    pipeWall.dimensions = Eigen::Vector3d {4,4,4};
    pipeWall.distance = 1;
    pipeWall.mass = -std::numeric_limits<double>::infinity();
    pipeWall.start_velocity = Eigen::Vector3d {0,0,0};

    Body liquid;
    liquid.shape = cuboid;
    liquid.fixpoint = Eigen::Vector3d{4,0,0};
    liquid.dimensions = Eigen::Vector3d {4,4,4};
    liquid.distance = 1;
    liquid.mass = 1;
    liquid.start_velocity = Eigen::Vector3d {3,-1,0};


    std::list<Particle> buf{};
    ParticleGenerator::generateCuboid(pipeWall, 0, buf, 3, 1, 1);
    ParticleGenerator::generateCuboid(liquid, 0, buf, 3, 1, 1);
    std::vector<Particle> bufVec(buf.begin(), buf.end());


    ParticleContainer pc(bufVec, {10,10,10}, 1, {}, true);

    //T = sum_particles(m*<v,v>)/(#dims*#particles)
    Thermostat thermostat(pc, 1, 2, 3, 100,0, true, ThermoMode::pipe);
    pc.forAllParticles([&](Particle& p){
        if(p.getM()>0){
            ASSERT_DOUBLE_EQ(std::sqrt(3*3 + 1*1), p.getV().norm());
        }else{
            ASSERT_DOUBLE_EQ(p.getV().norm(), 0);
        }
    });
    ASSERT_DOUBLE_EQ(thermostat.computeCurrentTemp(), (liquid.mass * 3 * 3)/(3));

    thermostat.notify();
    thermostat.notify();

    ASSERT_DOUBLE_EQ(thermostat.computeCurrentTemp(), 1);
    pc.forAllParticles([&](Particle& p){
        if(p.getM()>0){
            ASSERT_DOUBLE_EQ(-1, p.getV()[1]);
            ASSERT_DOUBLE_EQ((p.getV()+Eigen::Vector3d {0,1,0}).norm(), std::sqrt(3.0));
        }else{
            ASSERT_DOUBLE_EQ(p.getV().norm(), 0);
        }
    });
}