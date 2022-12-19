//
// Created by alex on 30.11.22.
//

#include "BoundsHandler.h"
#include "io/output/Logging.h"

namespace sim::physics::bounds {
    BoundsHandler::BoundsHandler(bound_t let, bound_t rit, bound_t tot,
                                 bound_t bot, bound_t frt, bound_t ret,
                                 force::ForceFunctorBase &ff, double st, double et, double dt, double eps, double sig,
                                 ParticleContainer &pc) :
            handleLeft(generateBound<side_t::left>(let, ff, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret)),
            handleRight(generateBound<side_t::right>(rit, ff, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret)),
            handleTop(generateBound<side_t::top>(tot, ff, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret)),
            handleBottom(generateBound<side_t::bottom>(bot, ff, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret)),
            handleFront(generateBound<side_t::front>(frt, ff, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret)),
            handleRear(generateBound<side_t::rear>(ret, ff, st, et, dt, eps, sig, pc, let, rit, bot, tot, frt, ret)),
            periodicActive(false), particleContainer(pc), forceFunctor(ff) {
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
        particleContainer.updateCells();

        //construct halo - only need to do 3 sides as rest is symmetric
        if(handleLeft->isPeriodic()) handleLeft->generateHalo();
        if(handleTop->isPeriodic()) handleTop->generateHalo();
        if(handleFront->isPeriodic()) handleFront->generateHalo();

        //calculate force - now need to handle opposite 3 sides, as these are where the hP are
        if(handleRight->isPeriodic()) handleRight->calcHaloForce();
        if(handleBottom->isPeriodic()) handleBottom->calcHaloForce();
        if(handleRear->isPeriodic()) handleRear->calcHaloForce();
    }
} // sim::physics::bounds

