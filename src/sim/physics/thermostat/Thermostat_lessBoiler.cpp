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
                               std::unordered_map<unsigned long, unsigned long> &id_to_index,
                               auto&){
            for(auto [_,a]: id_to_index){
                v[3*a] = beta * v[3*a];
                v[3*a+1] = beta * v[3*a+1];
                v[3*a+2] = beta * v[3*a+2];
            }
        });
    }
    else{//ThermoMode::pipe
        pc.runOnActiveData([&](std::vector<double> &force,
                               std::vector<double> &oldForce,
                               std::vector<double> &x,
                               std::vector<double> &v,
                               std::vector<double> &m,
                               std::vector<int> &type,
                               unsigned long count,
                               std::vector<double> &eps,
                               std::vector<double> &sig,
                               std::unordered_map<unsigned long, unsigned long> &id_to_index,
                               auto&){
            for(auto [_,a]: id_to_index){
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
                                    std::unordered_map<unsigned long, unsigned long> &id_to_index,
                                    auto&){
            for(auto [_,a]: id_to_index){
                sum += std::max(m[a], 0.) * (v[3*a]*v[3*a] + v[3*a+1]*v[3*a+1] + v[3*a+2]*v[3*a+2]);
            }
        });
        return sum/(dims*static_cast<double>(pc.activeSize()));
    }else{  //ThermoMode::pipe

        double sum{0};
        pc.runOnActiveData([this, &sum](std::vector<double> &force,
                                  std::vector<double> &oldForce,
                                  std::vector<double> &x,
                                  std::vector<double> &v,
                                  std::vector<double> &m,
                                  std::vector<int> &type,
                                  unsigned long count,
                                  std::vector<double> &eps,
                                  std::vector<double> &sig,
                                  std::unordered_map<unsigned long, unsigned long> &id_to_index,
                                  auto&){
            double meanYDirection{0};
            for(auto [_,a]: id_to_index){
                if(m[a] >= 0){
                    meanYDirection += v[3 * a + 1];
                }
            }
            meanYDirection = meanYDirection / static_cast<double>(numberFlowingParticles);

            for(auto [_,a]: id_to_index){
                sum += std::max(m[a], 0.) * (v[3*a]*v[3*a] + (v[3*a+1] - meanYDirection) * (v[3 * a + 1] - meanYDirection) + v[3 * a + 2] * v[3 * a + 2]);
            }
        });
        return sum/static_cast<double>(numberFlowingParticles);
    }
}
#endif
