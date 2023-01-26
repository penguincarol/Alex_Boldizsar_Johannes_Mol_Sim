#include "Particle.h"
#include "ParticleGenerator.h"
#include "Body.h"
#include "utils/MaxwellBoltzmannDistribution.h"

#include <iostream>

namespace ParticleGenerator {
    int bodyID = 1; /**every body gets generated with a unique body id */
    int particleID = 0; /**every particle gets generated with a unique particleID*/

	void generateCuboid(struct Body& body, double v_bolz, std::list<Particle>& buffer, int dims, double sigma, double epsilon){ //thermal friction hardcoded to 0.1, is that what we want to do?
        //Maybe it would be more efficient to concatenate two vectors instead of placing one particle after another in the ParticleContainer
        if(body.shape != cuboid){
            io::output::loggers::general->error("generateCuboid does not work for shapes that aren't Cuboids");
        }
        if(body.mass == -std::numeric_limits<double>::infinity()){v_bolz = 0;}  //catch pipeWalls
        int typeID = getNextBodyID();
        for (size_t x = 0; x < body.dimensions[0]; x++)
        {
            for (size_t y = 0; y < body.dimensions[1]; y++)
            {
                for (size_t z = 0; z < body.dimensions[2]; z++)
                {
                    Eigen::Vector3d pos = body.fixpoint + (body.distance * Eigen::Vector3d(x,y,z));
                    auto v_tmp = maxwellBoltzmannDistributedVelocity(v_bolz, dims);
                    Eigen::Vector3d v{ v_tmp[0], v_tmp[1], v_tmp[2] };
                    buffer.emplace_back(pos, (body.start_velocity +  v), body.mass, typeID, getNextParticleID());
                    buffer.back().setSigma(sigma);
                    buffer.back().setEpsilon(epsilon);
                }
            }
        }
    }

    /**
     * Helper methods to increase readability of generateMembrane.
     * This method is not meant to be used on its own
     * @return
     */
    static Eigen::Vector3d membrComputeOffset(int planeFlag, double x0, double x1){
        if(planeFlag == 0){
            return Eigen::Vector3d {0, x0, x1};
        }else if(planeFlag == 1){
            return Eigen::Vector3d {x0, 0, x1};
        }else{
            return Eigen::Vector3d {x0, x1, 0};
        }
    }

    void generateMembrane(struct Body& body, double v_bolz, std::list<Particle>& buffer, std::list<Membrane>& membranes, int dims, double sigma, double epsilon){
        if(body.shape != membrane){
            io::output::loggers::general->error("generateMembrane does not work on shapes that aren't Membranes");
        }
        if(body.mass == -std::numeric_limits<double>::infinity()){v_bolz = 0;}  //catch pipeWalls (should be irrelevant but you never know)

        //we initialize a 2 dimensional field with a 3 dimensional input where one of the dimensions is one
        //determine which dimension is one and which other 2 dimensions you can use for initialization
        int typeID = getNextBodyID() | 0x80000000;
        int planeFlag;  //0: y-z plane 1: x-z plane 2: x-y plane

        double x0;
        double x1;
        if(body.dimensions[0] == 1){
            x0 = body.dimensions[1];
            x1 = body.dimensions[2];
            planeFlag = 0;
        }else if(body.dimensions[1] == 1){
            x0 = body.dimensions[0];
            x1 = body.dimensions[2];
            planeFlag = 1;
        }else if(body.dimensions[2] == 1){
            x0 = body.dimensions[0];
            x1 = body.dimensions[1];
            planeFlag = 2;
        }else{
            io::output::loggers::general->error("Membrane needs to be 2 dimensional. Since the input specified a three dimensional Membrane it gets interpreted as a 2d Membrane in the x-y-plane");
            x0 = body.dimensions[1];
            x1 = body.dimensions[2];
            planeFlag = 0;
        }

        //std::vector<std::vector<unsigned long>> membrNodes(x0, std::vector<unsigned long>(0, x1));
        std::vector<std::vector<unsigned long>> membrNodes(x0, std::vector<unsigned long>{});
        for(double i = 0; i < x0; i++){
            for(double j = 0; j < x1; j++){
                Eigen::Vector3d pos = body.fixpoint + (body.distance * membrComputeOffset(planeFlag, i, j));
                auto v_tmp = maxwellBoltzmannDistributedVelocity(v_bolz, dims);

                int particleID = getNextParticleID();

                Eigen::Vector3d v { v_tmp[0], v_tmp[1], v_tmp[2] };
                buffer.emplace_back(pos, (body.start_velocity +  v), body.mass, typeID, particleID);
                buffer.back().setSigma(sigma);
                buffer.back().setEpsilon(epsilon);

                //membrNodes[i][j] = particleID; //set membraneNode to the index corresponding to the particle that is in this position in the "grid graph"
                membrNodes[i].emplace_back(particleID); //set membraneNode to the index corresponding to the particle that is in this position in the "grid graph"
            }
        }

        membranes.emplace_back(body.springStrength, body.desiredDistance, membrNodes, body.pullEndTime, std::array<double,3>{body.pullForce[0],body.pullForce[1],body.pullForce[2]}, body.pullIndices);

    }

    //this implementation is basically cutting a sphere out of a cuboid. that is obviously fine and asymptotically not worse than the best solution but still..
    //there might be a smarter way to do this
    void generateSphere(struct Body& body, double v_bolz, std::list<Particle>& buffer, int dims, double sigma, double epsilon){
        struct Body bodycopy(body); //we do configurate some parameters so better copy it.. we don't want weird side effects
        if(body.shape != sphere){
            io::output::loggers::general->error("generateSphere does not work for shapes that aren't Spheres");
        }
        if(body.mass == -std::numeric_limits<double>::infinity()){v_bolz = 0;}  //catch pipeWalls (should be irrelevant but you never know)

        //configuration stuff:
            //this should be the case in a well-formed input anyway:
        bodycopy.dimensions[1] = bodycopy.dimensions[0];
        bodycopy.dimensions[2] = bodycopy.dimensions[0];

        if(dims == 2){
            bodycopy.dimensions[2] = 0;
        }

        int typeID = getNextBodyID();
        for (size_t x = 0; x <= bodycopy.dimensions[0]; x++)
        {
            for (size_t y = 0; y <= bodycopy.dimensions[1]; y++)
            {
                for (size_t z = 0; z <= bodycopy.dimensions[2]; z++)
                {
                    if(x*x + y*y + z*z <= bodycopy.dimensions[0]*bodycopy.dimensions[0]){  //bodycopy.dimensions[0] is radius measured in particle-distances

                        //we only need to generate 1/8 of a sphere 8 times.. which is what we are doing here... it looks very ugly and seems very inefficient
                        /*std::vector<Eigen::Vector3d> pos{bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(x,y,z))};
                        if(x!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(-x,  y, z)));}
                        if(y!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d( x, -y, z)));}
                        if(z!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d( x,  y,-z)));}
                        if(x!=0 && y!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(-x,-y,z)));}
                        if(x!=0 && z!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(-x,y,-z)));}
                        if(y!=0 && z!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(x,-y,-z)));}
                        if(x!=0 && y!=0 && z!=0){{pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(-x,-y,-z)));}}*/

                        std::vector<Eigen::Vector3d> pos{bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(x,y,z))};
                        if(x!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(-(double)x,          y,         z)));}
                        if(y!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(         x, -(double)y,         z)));}
                        if(z!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(         x,          y,-(double)z)));}
                        if(x!=0 && y!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(-(double)x,-(double)y,         z)));}
                        if(x!=0 && z!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(-(double)x,         y,-(double)z)));}
                        if(y!=0 && z!=0){pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(         x,-(double)y,-(double)z)));}
                        if(x!=0 && y!=0 && z!=0){{pos.emplace_back(bodycopy.fixpoint + (bodycopy.distance * Eigen::Vector3d(-(double)x,-(double)y,-(double)z)));}}

                        for(auto p : pos){
                            auto v_tmp = maxwellBoltzmannDistributedVelocity(v_bolz, dims);
                            Eigen::Vector3d v { v_tmp[0], v_tmp[1], v_tmp[2] };
                            buffer.emplace_back(p, (body.start_velocity +  v), body.mass, typeID, getNextParticleID());
                            buffer.back().setSigma(sigma);
                            buffer.back().setEpsilon(epsilon);
                        }
                    }
                }
            }
        }
    }


    void generateParticle(Eigen::Vector3d& x, Eigen::Vector3d& v, double m, std::list<Particle>& buffer, double sigma, double epsilon){
        buffer.emplace_back(x, v, m, getNextBodyID(), getNextParticleID());
        buffer.back().setSigma(sigma);
        buffer.back().setEpsilon(epsilon);
    }

    int getNextBodyID() {
        return bodyID++;
    }

    int getNextParticleID(){
        return particleID++;
    }

}