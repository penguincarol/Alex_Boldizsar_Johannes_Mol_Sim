//
// Created by alex on 25.11.2022.
//

#pragma once


namespace io::input {
    /**
    * Defines names in the map for arguments in the input loader
    * */
    enum names {
        //IO
        outputFilePath,
        outputFileName,
        logLevel,
            // CHECKPOINT
            enableCheckpointing,
            simLastIteration,
        // SIM GENERAL
        startTime,
        endTime,
        delta_t,
        sigma,
        epsilon,
        brown,
        dimensions,
        // FORCE
        forceCalculation,
        enableGrav,
        enableMembrane,
        enableOMP,
        gGrav,
            // LINKED CELL
            enableLinkedCell,
            rCutoff,
            boundingBox_X0,
            boundingBox_X1,
            boundingBox_X2,
            boundCondFront,
            boundCondRear,
            boundCondLeft,
            boundCondRight,
            boundCondTop,
            boundCondBottom,
        // VELOCITY
        velocityCalculation,
        // POSITION
        positionCalculation,
        // BENCHMARK
        benchmark,
        benchmarkType,
        benchMaxBodySize,
        benchIterationCount,
        // THERMO
        enableThermo,
        thermoTInit,
        thermoNTerm,
        thermoTTarget,
        thermoDelta_t,
        //MISC
        names_count [[maybe_unused]]
    };
} // io::input