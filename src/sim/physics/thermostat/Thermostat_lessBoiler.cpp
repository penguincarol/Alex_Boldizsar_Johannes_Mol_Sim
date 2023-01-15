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
        pc.runOnActiveData([&](std::vector<double> &force,
                               std::vector<double> &oldForce,
                               std::vector<double> &x,
                               std::vector<double> &v,
                               std::vector<double> &m,
                               std::vector<int> &type,
                               unsigned long count,
                               std::vector<double> &eps,
                               std::vector<double> &sig,
                               std::vector<unsigned long> &activeParticles){
            for(auto a: activeParticles){
                v[3*a] = beta * v[3*a];
                v[3*a+1] = beta * v[3*a+1];
                v[3*a+2] = beta * v[3*a+2];
            }
        });
    }else{//ThermoMode::pipe
        pc.runOnActiveData([&](std::vector<double> &force,
                               std::vector<double> &oldForce,
                               std::vector<double> &x,
                               std::vector<double> &v,
                               std::vector<double> &m,
                               std::vector<int> &type,
                               unsigned long count,
                               std::vector<double> &eps,
                               std::vector<double> &sig,
                               std::vector<unsigned long> &activeParticles){
            for(auto a: activeParticles){
                v[3*a] = beta * v[3*a];
                v[3*a+2] = beta * v[3*a+2];
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
        //pc.forAllParticles([&sum](Particle& p ){sum += p.getM() * (p.getX().dot(p.getX()));});
        pc.runOnActiveData([&sum](std::vector<double> &force,
                                    std::vector<double> &oldForce,
                                    std::vector<double> &x,
                                    std::vector<double> &v,
                                    std::vector<double> &m,
                                    std::vector<int> &type,
                                    unsigned long count,
                                    std::vector<double> &eps,
                                    std::vector<double> &sig,
                                    std::vector<unsigned long> &activeParticles){
        for(auto a: activeParticles){
            sum += std::max(m[a], 0.) * (v[3*a]*v[3*a] + v[3*a+1]*v[3*a+1] + v[3*a+2]*v[3*a+2]);
        }
        });
        return sum/(dims*static_cast<double>(pc.activeSize()));
    }else{  //ThermoMode::pipe
        double meanYdirection{0};
        size_t numberFlowingParticles{0};   //TODO: there is A LOT of room for improvement here because we already know that!!
        pc.runOnActiveData([&meanYdirection, &numberFlowingParticles](std::vector<double> &force,
                                  std::vector<double> &oldForce,
                                  std::vector<double> &x,
                                  std::vector<double> &v,
                                  std::vector<double> &m,
                                  std::vector<int> &type,
                                  unsigned long count,
                                  std::vector<double> &eps,
                                  std::vector<double> &sig,
                                  std::vector<unsigned long> &activeParticles){
            for(auto a: activeParticles){
                if(m[a] >= 0){
                    meanYdirection += v[3*a+1];
                    numberFlowingParticles++;
                }
            }
            meanYdirection = meanYdirection / numberFlowingParticles;   //unfortunately we can't use activeParticles.size() here because of pipe particles
        });

        double sum{0};
        //pc.forAllParticles([&sum](Particle& p ){sum += p.getM() * (p.getX().dot(p.getX()));});
        pc.runOnActiveData([&sum, meanYdirection](std::vector<double> &force,
                                  std::vector<double> &oldForce,
                                  std::vector<double> &x,
                                  std::vector<double> &v,
                                  std::vector<double> &m,
                                  std::vector<int> &type,
                                  unsigned long count,
                                  std::vector<double> &eps,
                                  std::vector<double> &sig,
                                  std::vector<unsigned long> &activeParticles){
            for(auto a: activeParticles){
                sum += std::max(m[a], 0.) * (v[3*a]*v[3*a] + (v[3*a+1]-meanYdirection)*(v[3*a+1]-meanYdirection) + v[3*a+2]*v[3*a+2]);
            }
        });
        return sum/numberFlowingParticles;

    }
}
#endif
