//
// Created by johnny on 09.01.23.
//

#include "FMembrane.h"

namespace sim::physics::force {

    /*
    void FMembrane::operator()() {
        particleContainer.forAllMembraneSprings([&](Particle &p1, Particle &p2, double desiredDistance, double springStrength){
            //this method is basically the pairFun but you can't really make it pretty with the signature of pairFun..
           Eigen::Vector3d dist{p2.getX() - p1.getX()};   //p1 = pi, p2 = pj
           double norm = dist.norm();
           Eigen::Vector3d f_ij = springStrength * (norm - desiredDistance) * dist/norm;
           p1.add_to_F(f_ij);
           p2.add_to_F(-f_ij);
        });
    }*/

    void FMembrane::addSpringForce(size_t p1i, size_t p1j, size_t p2i, size_t p2j,
                                           Membrane& membrane, std::vector<double>& force, std::vector<double>& x){
        size_t idI = membrane.getMembrNodes()[p1i][p1j];
        size_t idJ = membrane.getMembrNodes()[p2i][p2j];

        //if that's not how ids work i have to rewrite this part but that should be rather easy
        //TODO: x[pc.idMap[idI]*3 + 0] instead of what is currently there

        Eigen::Vector3d dist{x[idI*3 + 0] - x[idJ*3 + 0],
                             x[idI*3 + 1] - x[idJ*3 + 1],
                             x[idI*3 + 2] - x[idJ*3 + 2]};
        double norm = dist.norm();
        Eigen::Vector3d f_ij = membrane.getSpringStrength() * (norm - membrane.getDesiredDistance()) * dist / norm;

        force[idI*3 + 0] -= f_ij[0];
        force[idI*3 + 1] -= f_ij[1];
        force[idI*3 + 2] -= f_ij[2];

        force[idJ*3 + 0] += f_ij[0];
        force[idJ*3 + 1] += f_ij[1];
        force[idJ*3 + 2] += f_ij[2];
    }

    void FMembrane::operator()() {
        //fun(membranes, force, x, count);
        particleContainer.runOnMembranes([&](std::vector<Membrane>& membranes,
                                        std::vector<double>& force,
                                        std::vector<double>& x,
                                        unsigned long count){
            for(Membrane& membrane:membranes){
                std::vector<std::vector<unsigned long>>& grid = membrane.getMembrNodes();
                if(grid.size() <= 0){continue;}


                //actually trying to avoid code duplication this time
                const size_t num_cases{4};
                std::array<size_t, num_cases> lowerboundsI = {0,0,0,0};
                std::array<size_t, num_cases> lowerboundsJ = {0,0,0,1};

                std::array<size_t, num_cases> upperboundsI = {grid.size()-1, grid.size(), grid.size()-1, grid.size()-1};
                std::array<size_t, num_cases> upperboundsJ = {grid[0].size(), grid[0].size()-1, grid[0].size()-1, grid[0].size()};

                using a = std::array<int, 2>;
                std::array<a, num_cases> offsets = {a{1,0}, //horizontal
                                            a{0,1}, //vertical
                                            a{1,1}, //top right
                                            a{1,-1} //bottom right
                                            };

                for(size_t _case = 0; _case < num_cases; _case++){

                    for(size_t i = lowerboundsI[_case]; i < upperboundsI[_case]; i++){
                        for(size_t j = lowerboundsJ[_case]; j < upperboundsJ[_case]; j++){
                            this->addSpringForce(i,j, i+ offsets[_case][0], j + offsets[_case][1],
                                                                      membrane, force, x);

                        }
                    }

                }

            }
        });
    }

    void FMembrane::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
    }

    pair_fun_t &FMembrane::getForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};;
    }

    fpair_fun_t FMembrane::getFastForceFunction() {
        throw std::runtime_error{"This should not be called. Not supported."};
    }
} // sim::physics::force