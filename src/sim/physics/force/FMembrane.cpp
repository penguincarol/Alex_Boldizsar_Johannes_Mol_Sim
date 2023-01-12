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

    void FMembrane::operator()() {
        forceDelegate->operator()();    //take it out for new forceFunctor-interpretation

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
                            //TODO: computeSpringForces
                        }
                    }

                }

                //this part is getting redundant because of the no-duplication stuff going on above
                //horizontal lines
                for(size_t i = 0; i < grid.size()-1; i++){
                    for(size_t j = 0; j < grid[0].size(); j++){
                        //TODO: computeSpringForces
                    }
                }

                //vertical lines
                for(size_t i = 0; i < grid.size(); i++){
                    for(size_t j = 0; j < grid[0].size()-1; j++){
                        //TODO: computeSpringForces
                    }
                }

                //to top right diagonals
                for(size_t i = 0; i < grid.size()-1; i++){
                    for(size_t j = 0; j < grid[0].size()-1; j++){
                        //TODO: computeSpringForces
                    }
                }

                //to bottom right diagonals
                for(size_t i = 0; i < grid.size()-1; i++){
                    for(size_t j = 1; j < grid[0].size(); j++){
                        //TODO: computeSpringForces
                    }
                }
                //end of redundant stuff


            }
        });
    }

    pair_fun_t &FMembrane::getForceFunction() {
        return pairFun;
    }

    static void fastPairFunction(std::vector<double> &force,
                                 std::vector<double> &x,
                                 std::vector<double> &eps,
                                 std::vector<double> &sig,
                                 std::vector<double> &m,
                                 unsigned long indexI, unsigned long indexJ, bool wbI, bool wbJ) {
        //TODO (if i find a way to make it pretty or we need to change signatures)
    }

    void FMembrane::setPairFun() {
        pairFun = [&](Particle &p1, Particle &p2) {
            //TODO (if i find a way to make it pretty or we need to change signatures)
        };
        fpairFun = fastPairFunction;
    }

    void FMembrane::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
        forceDelegate->setParticleContainer(pc);
        setPairFun();
        membranes = pc.getMembranes();
    }

    fpair_fun_t FMembrane::getFastForceFunction() {
        return fpairFun;
    }
} // sim::physics::force