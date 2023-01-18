#include "Thermostat.h"
#include "data/ParticleContainer.h"
#include "data/Particle.h"
#include "utils/MaxwellBoltzmannDistribution.h"

#include <vector>
#include <cmath>

void Thermostat::setParticleContainer(ParticleContainer& particleContainer){
    this->pc = particleContainer;
}

double Thermostat::computeBeta(){
    double Tcur{computeCurrentTemp()};

    double Tnew;
    if((Ttarget - Tcur) > deltaTemp){
        Tnew = Tcur + deltaTemp;
    }else if((Ttarget-Tcur) < -deltaTemp){
        Tnew = Tcur - deltaTemp;
    }else{
        Tnew = this->Ttarget;
    }

    return std::sqrt(Tnew/Tcur);
}

void Thermostat::notify(){
    //TODO: depending on the answer (see discussion in discord) isActive needs to get deleted as a variable 
    //and the if(isActive)-part needs to be taken into the first if-statement
    //or this part is correct so far..
    countSinceLastActivation++;
    if(countSinceLastActivation >= countThreshold){
        countSinceLastActivation -= countThreshold;

        getCooking();
    }
}


void Thermostat::initializeBrownTemp(double TInit){
    pc.runOnData([&](std::vector<Particle>& particles,
                                    std::vector<Membrane>& membranes,
                                    ParticleContainer::VectorCoordWrapper& cells,
                                    unsigned long count,
                                    std::vector<unsigned long>& activeParticles,
                                    std::unordered_map<unsigned long, unsigned long> &id_to_index){
        for(auto [_,a]: id_to_index){
            auto brown{maxwellBoltzmannDistributedVelocity(std::sqrt(TInit/particles[a].getM()), dims)};
            //std::array<double, 3> brown{0.,0.,0.};
            particles[a].add_to_V({brown[0],brown[1],brown[2]});
        }
    });
}