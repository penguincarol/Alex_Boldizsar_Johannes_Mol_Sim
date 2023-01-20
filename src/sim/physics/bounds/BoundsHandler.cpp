//
// Created by alex on 30.11.22.
//

#include "BoundsHandler.h"
#include "io/output/Logging.h"

namespace sim::physics::bounds {
    BoundsHandler::BoundsHandler(bound_t let, bound_t rit, bound_t tot,
                                 bound_t bot, bound_t frt, bound_t ret,
                                 force::ForceHandler &fh, double st, double et, double dt, double eps, double sig,
                                 ParticleContainer &pc, bool eOMP) :
            handleLeft(generateBound<side_t::left>(let, fh, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret, eOMP)),
            handleRight(generateBound<side_t::right>(rit, fh, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret, eOMP)),
            handleTop(generateBound<side_t::top>(tot, fh, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret, eOMP)),
            handleBottom(generateBound<side_t::bottom>(bot, fh, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret, eOMP)),
            handleFront(generateBound<side_t::front>(frt, fh, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret, eOMP)),
            handleRear(generateBound<side_t::rear>(ret, fh, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret, eOMP)),
            periodicActive(false), particleContainer(pc), enableOMP(eOMP) {
        //check for null pointers
        if (handleLeft == nullptr || handleRight == nullptr || handleTop == nullptr ||
            handleBottom == nullptr || handleFront == nullptr || handleRear == nullptr) {
            io::output::loggers::general->error("Failed to create Bounds Handler! Malloc Failed.");
            exit(-1);
        }

        //check for periodic
        if(let == periodic || rit == periodic || tot == periodic || bot == periodic || frt == periodic || ret == periodic)
            periodicActive = true;
    }

    BoundsHandler::~BoundsHandler() {
        delete handleLeft;
        delete handleRight;
        delete handleTop;
        delete handleBottom;
        delete handleFront;
        delete handleRear;
    }

    void BoundsHandler::setParticleContainer(ParticleContainer &pc) {
        handleLeft->setParticleContainer(pc);
        handleRight->setParticleContainer(pc);
        handleTop->setParticleContainer(pc);
        handleBottom->setParticleContainer(pc);
        handleFront->setParticleContainer(pc);
        handleRear->setParticleContainer(pc);
    }

    void BoundsHandler::operator()() {
        if(periodicActive) handlePeriodic();

        if(!handleLeft->isPeriodic()) handleLeft->operator()();
        if(!handleRight->isPeriodic()) handleRight->operator()();
        if(!handleTop->isPeriodic()) handleTop->operator()();
        if(!handleBottom->isPeriodic()) handleBottom->operator()();
        if(!handleFront->isPeriodic()) handleFront->operator()();
        if(!handleRear->isPeriodic()) handleRear->operator()();
    }

    void BoundsHandler::handlePeriodic() {
        particleContainer.clearHalo();
        //move on border
        if(handleLeft->isPeriodic()) handleLeft->operator()();
        if(handleRight->isPeriodic()) handleRight->operator()();
        if(handleTop->isPeriodic()) handleTop->operator()();
        if(handleBottom->isPeriodic()) handleBottom->operator()();
        if(handleFront->isPeriodic()) handleFront->operator()();
        if(handleRear->isPeriodic()) handleRear->operator()();
        //particleContainer.updateCells();

        //construct halo - only need to do 3 sides as rest is symmetric
        if(handleLeft->isPeriodic()) handleLeft->generateHalo();
        if(handleBottom->isPeriodic()) handleBottom->generateHalo();
        if(handleFront->isPeriodic()) handleFront->generateHalo();

        //calculate force - now need to handle opposite 3 sides, as these are where the hP are
        if(handleRight->isPeriodic()) handleRight->calcHaloForce();
        if(handleTop->isPeriodic()) handleTop->calcHaloForce();
        if(handleRear->isPeriodic()) handleRear->calcHaloForce();
    }
} // sim::physics::bounds

