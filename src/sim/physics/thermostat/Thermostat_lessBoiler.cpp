#include "Thermostat.h"
#include "data/ParticleContainer.h"
#include "data/Particle.h"
#include "io/output/Logging.h"

#include <vector>


#ifndef slow
void Thermostat::getCooking(){
    io::output::loggers::simulation->debug("getCooking called");
    double beta{computeBeta()};

    //v = beta*v for all active Particles
    if(thermoMode == ThermoMode::normal){
        pc.runOnData([&](std::vector<Particle>& particles,
                                        std::vector<Membrane>& membranes,
                                        ParticleContainer::VectorCoordWrapper& cells,
                                        unsigned long count,
                                        std::vector<unsigned long>& activeParticles,
                                        std::unordered_map<unsigned long, unsigned long> &id_to_index){
            for(auto [_,a]: id_to_index){
                particles[a].setV(particles[a].getV() * beta);
            }
        });
    }
    else{//ThermoMode::pipe
        pc.runOnData([&](std::vector<Particle>& particles,
                         std::vector<Membrane>& membranes,
                         ParticleContainer::VectorCoordWrapper& cells,
                         unsigned long count,
                         std::vector<unsigned long>& activeParticles,
                         std::unordered_map<unsigned long, unsigned long> &id_to_index){
            for(auto [_,a]: id_to_index){
                particles[a].setV({
                    beta * particles[a].getV()[0],
                           particles[a].getV()[1],
                    beta * particles[a].getV()[2],
                });
            }
        });
    }
    io::output::loggers::simulation->debug("The temperature after letting the thermostat work is " + std::to_string(computeCurrentTemp()));
}


double Thermostat::computeCurrentTemp(){
    //E_kin = sum_particles (m* <v,v>/2)
    //E_kin = #dims*#particles*T/2

    //T = sum_particles(m*<v,v>)/(#dims*#particles)

    if(thermoMode == ThermoMode::normal){
        double sum{0};
        pc.forAllParticles([&sum](Particle& p ){sum += p.getM() * (p.getX().dot(p.getX()));});
        return sum/(dims*static_cast<double>(pc.activeSize()));
    }else{  //ThermoMode::pipe

        double sum{0};
        pc.runOnData([&](std::vector<Particle>& particles,
                         std::vector<Membrane>& membranes,
                         ParticleContainer::VectorCoordWrapper& cells,
                         unsigned long count,
                         std::vector<unsigned long>& activeParticles,
                         std::unordered_map<unsigned long, unsigned long> &id_to_index){
            double meanYDirection{0};
            for(auto [_,a]: id_to_index){
                if(particles[a].getM() >= 0){
                    meanYDirection += particles[a].getV()[1];
                }
            }
            meanYDirection = meanYDirection / static_cast<double>(numberFlowingParticles);

            for(auto [_,a]: id_to_index){
                Particle& p = particles[a];
                sum += std::max(p.getM(), 0.) *
                        (p.getV()[0]*p.getV()[0] + (p.getV()[1] - meanYDirection) * (p.getV()[1] - meanYDirection) + p.getV()[2] * p.getV()[2]);
            }
        });

        return sum/static_cast<double>(numberFlowingParticles);
    }
}
#endif
