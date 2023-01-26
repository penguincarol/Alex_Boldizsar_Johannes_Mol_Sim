//
// Created by johnny on 09.01.23.
//

#include "FMembrane.h"

namespace sim::physics::force {

    void FMembrane::addSpringForce(size_t p1i, size_t p1j, size_t p2i, size_t p2j,
                                           Membrane& membrane, Kokkos::View<double*>& force, Kokkos::View<double*>& x, std::unordered_map<unsigned long, unsigned long> &id_to_index){
        size_t idI = membrane.getMembrNodes()[p1i][p1j];
        size_t idJ = membrane.getMembrNodes()[p2i][p2j];
        size_t indexI = id_to_index[idI];
        size_t indexJ = id_to_index[idJ];

        double d0,d1,d2;
        d0 = x[indexI * 3 + 0] - x[indexJ * 3 + 0];
        d1 = x[indexI * 3 + 1] - x[indexJ * 3 + 1];
        d2 = x[indexI * 3 + 2] - x[indexJ * 3 + 2];
        double f0,f1,f2;
        double norm = std::sqrt(d0*d0+d1*d1+d2*d2);
        double k = membrane.getSpringStrength();
        double dd = membrane.getDesiredDistance();
        if(p1i != p2i && p1j != p2j){
            f0 = k * (norm - std::sqrt(2.0) * dd) * d0 / norm;
            f1 = k * (norm - std::sqrt(2.0) * dd) * d1 / norm;
            f2 = k * (norm - std::sqrt(2.0) * dd) * d2 / norm;
        }else{
            f0 = k * (norm - dd) * d0 / norm;
            f1 = k * (norm - dd) * d1 / norm;
            f2 = k * (norm - dd) * d2 / norm;
        }

        force[indexI * 3 + 0] -= f0;
        force[indexI * 3 + 1] -= f1;
        force[indexI * 3 + 2] -= f2;

        force[indexJ * 3 + 0] += f0;
        force[indexJ * 3 + 1] += f1;
        force[indexJ * 3 + 2] += f2;
    }

    void FMembrane::operator()() {
        //fun(membranes, force, x, count);
        particleContainer.runOnMembranes([](std::vector<Membrane>& membranes,
                                            Kokkos::View<double*>& force,
                                            Kokkos::View<double*>& x,
                                            unsigned long count,
                                            std::unordered_map<unsigned long, unsigned long> &id_to_index){
            for(Membrane& membrane:membranes){
                auto& grid = membrane.getMembrNodes();
                if(grid.empty()) continue;

                for(size_t i = 0; i < grid.size()-1; i++) {
                    for(size_t j = 0; j < grid[i].size(); j++)
                        addSpringForce(i,j,i+1,j,membrane,force,x,id_to_index);
                }
                for(size_t i = 0; i < grid.size(); i++) {
                    for(size_t j = 0; j < grid[i].size()-1; j++)
                        addSpringForce(i,j,i,j+1,membrane,force,x,id_to_index);
                }
                for(size_t i = 0; i < grid.size()-1; i++) {
                    for(size_t j = 0; j < grid[i].size()-1; j++)
                        addSpringForce(i,j,i+1,j+1,membrane,force,x,id_to_index);
                }
                for(size_t i = 0; i < grid.size()-1; i++) {
                    for(size_t j = 1; j < grid[i].size(); j++)
                        addSpringForce(i,j,i+1,j-1,membrane,force,x,id_to_index);
                }
            }
        });
    }

    void FMembrane::setParticleContainer(ParticleContainer &pc) {
        PhysicsFunctorBase::setParticleContainer(pc);
    }

} // sim::physics::force