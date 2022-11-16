#include "ParticleContainer.h"
#include "Particle.h"
#include "io/Logging.h"

#include <vector>
#include <thread>
#include <Eigen>

ParticleContainer::ParticleContainer() = default;

ParticleContainer::ParticleContainer(const std::vector<Particle>& buffer) {
    for (const auto& p : buffer) particles.emplace_back(p);
}

unsigned long ParticleContainer::size() {
    return particles.size();
}

std::vector<Particle> &ParticleContainer::getParticles() {
    return particles;
}

Particle &ParticleContainer::getParticle(unsigned long i) {
    if (i < particles.size()) {
        return particles[i];
    } else {
        throw std::runtime_error("Tried to access Particle with index out of bounds!\n");
    }
}

void ParticleContainer::forAllParticles(void (function)(Particle &p)) {
    std::for_each(particles.begin(), particles.end(), [&](Particle &g) { function(g); });
}

void ParticleContainer::forAllPairs(void (function)(Particle &p1, Particle &p2)) {
    for (u_int32_t i = 0; i < particles.size(); i++) {
        for (u_int32_t j = i + 1; j < particles.size(); j++) {
            Particle &p1 = particles[i];
            Particle &p2 = particles[j];
            function(p1, p2);
        }
    }
}

void ParticleContainer::forceInIndexSpan(size_t start, size_t end, std::vector<Eigen::Vector3d>& acc, 
            Eigen::Vector3d (function)(Particle &p1, Particle &p2), std::vector<Particle> & particles){
    for (u_int32_t i = start; i < end; i++) {
        for (u_int32_t j = i + 1; j < particles.size(); j++) {
            Particle &p1 = particles[i];
            Particle &p2 = particles[j];
            Eigen::Vector3d force{function(p1, p2)};
            acc[i] += force;
            acc[j] += -force;
        }
    }
}

void ParticleContainer::addUpForces(Eigen::Vector3d (function)(Particle &p1, Particle &p2)){
    int number_of_threads = std::thread::hardware_concurrency();
    loggers::general->debug("Creating {} threads to add up Forces", number_of_threads);
    std::vector<Eigen::Vector3d> force_acc(particles.size(), {0.,0.,0.});
    std::vector<std::vector<Eigen::Vector3d>> force_accs;
    std::vector<std::thread> threads;

    //define ranges that are to be distributed;
    std::vector<size_t> splits;
    for(int i = 0; i < number_of_threads; i++){splits.emplace_back(i* particles.size() / number_of_threads);};
    splits.emplace_back(particles.size());  //we need n+1 splitpoints to define n ranges in between that are to be distributed

    for(int i = 0; i < number_of_threads; i++){
        threads.emplace_back(std::thread(forceInIndexSpan ,splits[i], splits[i+1], std::ref(force_accs[i]), std::ref(function), std::ref(particles)));};

    for(std::thread & t : threads){
        if (t.joinable())
            t.join();
    }

    for(std::vector<Eigen::Vector3d> & cur_acc : force_accs){
        for(size_t i = 0; i < particles.size(); i++){
            particles[i].add_to_F(cur_acc[i]);
        }
    }
}

ParticleContainer::Iterator ParticleContainer::begin() {
    return Iterator{particles.begin()} ;
}

ParticleContainer::Iterator ParticleContainer::end() {
    return Iterator{particles.end()};
}

void ParticleContainer::addParticle(const Particle& p){
    particles.emplace_back(p);
}

void ParticleContainer::clear(){
    particles.clear();
}
