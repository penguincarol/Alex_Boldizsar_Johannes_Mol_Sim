//
// Created by alex on 30.11.22.
//

#pragma once

#include <string>
#include <variant>
#include <iostream>
#include "sim/physics/bounds/types.h"
#include "sim/physics/force/types.h"
#include "sim/physics/velocity/types.h"
#include "sim/physics/position/types.h"
#include "arg_names.h"
#include "io/input/cli/CLIArgs.h"
#include "sim/physics/thermostat/Thermostat.h"

namespace io::input {

    using force_t = sim::physics::force::type;
    using pos_t = sim::physics::position::type;
    using vel_t = sim::physics::velocity::type;
    using bound_t = sim::physics::bounds::type;

    struct ERROR_T {
    };


    /**
    * Defines a map from enum name to new data type
    * */
    template<io::input::names N>
    struct map {
        struct ERROR_T {
        };
        struct UNKNOWN_T {
            using key = ERROR_T;
            using value = ERROR_T;
        };
        struct WRONG_ORDER_T {
            using key = ERROR_T;
            using value = ERROR_T;
        };
        template<io::input::names K, typename V>
        struct entry {
            static constexpr io::input::names key = K;
            using value = V;
        };

        template<bool b, typename T, typename F>
        struct cond {
            using type = T;
        };
        template<typename T, typename F>
        struct cond<false, T, F> {
            using type = F;
        };

        template<int count, typename E1 = UNKNOWN_T, typename... E>
        struct extract {
            using type = typename cond<E1::key==N, E1, typename extract<count-1, E...>::type>::type;
        };

        template<typename E1, typename... E>
        struct extract<1, E1, E...>{
            using type = typename extract<0, E...>::type;
        };

        template<typename E1, typename... E>
        struct extract<0, E1, E...> {
            using type = E1;
        };

        //Lookup table
        //Order of entries must be the same as in io::input::names
        using current_entry = typename extract<static_cast<int>(N),
                entry<outputFilePath, std::string>,
                entry<outputFileName, std::string>,
                entry<logLevel, int>,
                entry<enableCheckpointing, bool>,
                entry<simLastIteration, int>,
                entry<startTime, double>,
                entry<endTime, double>,
                entry<delta_t, double>,
                entry<sigma, double>,
                entry<epsilon, double>,
                entry<brown, double>,
                entry<dimensions, int>,
                entry<forceCalculation, force_t>,
                entry<enableGrav, bool>,
                entry<enableMembrane, bool>,
                entry<enableMembranePull, bool>,
                entry<enableOMP, bool>,
                entry<enableProfiler, bool>,
                entry<profilerNumBins, int>,
                entry<gGrav0, double>,
                entry<gGrav1, double>,
                entry<gGrav2, double>,
                entry<enableLinkedCell, bool>,
                entry<rCutoff, double>,
                entry<boundingBox_X0, double>,
                entry<boundingBox_X1, double>,
                entry<boundingBox_X2, double>,
                entry<boundCondFront, bound_t>,
                entry<boundCondRear, bound_t>,
                entry<boundCondLeft, bound_t>,
                entry<boundCondRight, bound_t>,
                entry<boundCondTop, bound_t>,
                entry<boundCondBottom, bound_t>,
                entry<velocityCalculation, vel_t>,
                entry<positionCalculation, pos_t>,
                entry<benchmark, bool>,
                entry<benchmarkType, std::string>,
                entry<benchMaxBodySize, int>,
                entry<benchIterationCount, int>,
                entry<enableThermo, bool>,
                entry<thermoTInit, double>,
                entry<thermoNTerm, int>,
                entry<thermoTTarget, double>,
                entry<thermoDelta_t, double>,
                entry<thermoType_t, ThermoMode>
        >::type;

        using get = typename cond<current_entry::key != N, WRONG_ORDER_T, typename current_entry::value>::type;
    };

    /**
     * Manages all input parameters. Will first acquire params from CLI and then File.
     * */
    class Configuration {

    private:
        //value lock
        std::unordered_map<names, bool> valueLock;

        //data storage
        using config_entry_t = std::variant<std::string, double, int, bool, force_t, pos_t, vel_t, bound_t, ThermoMode, ERROR_T>;
        std::unordered_map<names, config_entry_t> dataStorage;

    public:

        /**
         * Gets the value of field N, which is of io::input::names, and has the corresponding data type.
         * */
        template<names N, typename R = typename map<N>::get>
        R get() {
            return std::get<R>(dataStorage[N]);
        }

        /**
         * Returns the internal data map. For testing purposes only
         * */
        std::unordered_map<names, config_entry_t> &getData() {
            return dataStorage;
        }

        /**
         * Returns the internal lock map. For testing purposes only
         * */
        std::unordered_map<names, bool> &getLocks() {
            return valueLock;
        }

        /**
         * Loads in all args from the cli registry.
         * CLIArgsParser::parseArgs() has to be called before.
         * */
        void loadCLIArgs() {
            dataStorage[outputFilePath] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-of")).value;
            dataStorage[outputFileName] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-o")).value;
            dataStorage[startTime] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-st")).value;
            dataStorage[endTime] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-et")).value;
            dataStorage[delta_t] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-dt")).value;
            dataStorage[forceCalculation] = sim::physics::force::stot(std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-f")).value);
            dataStorage[positionCalculation] = sim::physics::position::stot(std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-x")).value);
            dataStorage[velocityCalculation] = sim::physics::velocity::stot(std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-v")).value);
            dataStorage[sigma] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-sig")).value;
            dataStorage[epsilon] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-eps")).value;
            dataStorage[brown] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-brown")).value;
            dataStorage[enableLinkedCell] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-fELC")).value != 0;
            dataStorage[enableGrav] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-fEGrav")).value != 0;
            dataStorage[enableMembrane] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-fEMem")).value != 0;
            dataStorage[enableMembranePull] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-fEMemPull")).value != 0;
            dataStorage[enableOMP] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-fEOMP")).value != 0;
            dataStorage[rCutoff] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-rc")).value;
            dataStorage[boundingBox_X0] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-bbox0")).value;
            dataStorage[boundingBox_X1] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-bbox1")).value;
            dataStorage[boundingBox_X2] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-bbox2")).value;
            dataStorage[boundCondFront] = sim::physics::bounds::stot(std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndFront")).value);
            dataStorage[boundCondRear] = sim::physics::bounds::stot(std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndRear")).value);
            dataStorage[boundCondLeft] = sim::physics::bounds::stot(std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndLeft")).value);
            dataStorage[boundCondRight] = sim::physics::bounds::stot(std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndRight")).value);
            dataStorage[boundCondTop] = sim::physics::bounds::stot(std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndTop")).value);
            dataStorage[boundCondBottom] = sim::physics::bounds::stot(std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndBottom")).value);
            dataStorage[dimensions] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-dims")).value;
            dataStorage[logLevel] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-llv")).value;
            dataStorage[benchmark] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bench")).isSet;
            dataStorage[benchmarkType] = std::get<io::input::ArgEntry<std::string>>(io::input::cli_arg_map.at("-bench")).value;
            dataStorage[benchMaxBodySize] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-bMax")).value;
            dataStorage[benchIterationCount] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-i")).value;
            dataStorage[enableThermo] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-thermo")).value != 0;
            dataStorage[thermoTInit] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-ti")).value;
            dataStorage[thermoNTerm] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-nt")).value;
            dataStorage[thermoTTarget] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-tt")).value;
            dataStorage[thermoDelta_t] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-dTemp")).value;
            dataStorage[enableCheckpointing] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-cp")).value != 0;
            dataStorage[gGrav0] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-gGrav0")).value;
            dataStorage[gGrav1] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-gGrav1")).value;
            dataStorage[gGrav2] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-gGrav2")).value;
            dataStorage[simLastIteration] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-lastIt")).value;
            dataStorage[enableProfiler] = false;
            dataStorage[profilerNumBins] = 0;
            dataStorage[thermoType_t] = ThermoMode::normalMode;

            valueLock[outputFilePath] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-of")).isSet;
            valueLock[outputFileName] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-o")).isSet;
            valueLock[startTime] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-st")).isSet;
            valueLock[endTime] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-et")).isSet;
            valueLock[delta_t] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-dt")).isSet;
            valueLock[forceCalculation] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-f")).isSet;
            valueLock[positionCalculation] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-x")).isSet;
            valueLock[velocityCalculation] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-v")).isSet;
            valueLock[sigma] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-sig")).isSet;
            valueLock[epsilon] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-eps")).isSet;
            valueLock[brown] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-brown")).isSet;
            valueLock[enableLinkedCell] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-fELC")).isSet;
            valueLock[enableGrav] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-fEGrav")).isSet;
            valueLock[enableMembrane] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-fEMem")).isSet;
            valueLock[enableMembranePull] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-fEMemPull")).isSet;
            valueLock[enableOMP] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-fEOMP")).isSet;
            valueLock[rCutoff] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-rc")).isSet;
            valueLock[boundingBox_X0] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-bbox0")).isSet;
            valueLock[boundingBox_X1] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-bbox1")).isSet;
            valueLock[boundingBox_X2] = std::get<io::input::ArgEntry<double>>(cli_arg_map.at("-bbox2")).isSet;
            valueLock[boundCondFront] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndFront")).isSet;
            valueLock[boundCondRear] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndRear")).isSet;
            valueLock[boundCondLeft] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndLeft")).isSet;
            valueLock[boundCondRight] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndRight")).isSet;
            valueLock[boundCondTop] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndTop")).isSet;
            valueLock[boundCondBottom] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bndBottom")).isSet;
            valueLock[dimensions] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-dims")).isSet;
            valueLock[logLevel] = std::get<io::input::ArgEntry<int>>(cli_arg_map.at("-llv")).isSet;
            valueLock[benchmark] = std::get<io::input::ArgEntry<std::string>>(cli_arg_map.at("-bench")).isSet;
            valueLock[benchmarkType] = std::get<io::input::ArgEntry<std::string>>(io::input::cli_arg_map.at("-bench")).isSet;
            valueLock[benchMaxBodySize] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-bMax")).isSet;
            valueLock[benchIterationCount] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-i")).isSet;
            valueLock[enableThermo] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-thermo")).isSet;
            valueLock[thermoTInit] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-ti")).isSet;
            valueLock[thermoNTerm] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-nt")).isSet;
            valueLock[thermoTTarget] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-tt")).isSet;
            valueLock[thermoDelta_t] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-dTemp")).isSet;
            valueLock[enableCheckpointing] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-cp")).isSet;
            valueLock[gGrav0] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-gGrav0")).isSet;
            valueLock[gGrav1] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-gGrav1")).isSet;
            valueLock[gGrav2] = std::get<io::input::ArgEntry<double>>(io::input::cli_arg_map.at("-gGrav2")).isSet;
            valueLock[simLastIteration] = std::get<io::input::ArgEntry<int>>(io::input::cli_arg_map.at("-lastIt")).isSet;
            valueLock[enableProfiler] = false;
            valueLock[profilerNumBins] = false;
            valueLock[thermoType_t] = false;
        }

        /**
         * Loads in all args from the provided arg map. This should be from the IOWrapper.
         * */
        void loadIOWArgs(const std::unordered_map<names, std::string> &argMap) {
            if (!valueLock[outputFilePath] && argMap.contains(outputFilePath))dataStorage[outputFilePath] = argMap.at(outputFilePath);
            if (!valueLock[outputFileName] && argMap.contains(outputFileName))dataStorage[outputFileName] = argMap.at(outputFileName);
            if (!valueLock[startTime] && argMap.contains(startTime))dataStorage[startTime] = std::stod(argMap.at(startTime));
            if (!valueLock[endTime] && argMap.contains(endTime)) dataStorage[endTime] = std::stod(argMap.at(endTime));
            if (!valueLock[delta_t] && argMap.contains(delta_t)) dataStorage[delta_t] = std::stod(argMap.at(delta_t));
            if (!valueLock[forceCalculation] && argMap.contains(forceCalculation))dataStorage[forceCalculation] = sim::physics::force::stot(argMap.at(forceCalculation));
            if (!valueLock[positionCalculation] && argMap.contains(positionCalculation))dataStorage[positionCalculation] = sim::physics::position::stot(argMap.at(positionCalculation));
            if (!valueLock[velocityCalculation] && argMap.contains(velocityCalculation))dataStorage[velocityCalculation] = sim::physics::velocity::stot(argMap.at(velocityCalculation));
            if (!valueLock[sigma] && argMap.contains(sigma)) dataStorage[sigma] = std::stod(argMap.at(sigma));
            if (!valueLock[epsilon] && argMap.contains(epsilon)) dataStorage[epsilon] = std::stod(argMap.at(epsilon));
            if (!valueLock[brown] && argMap.contains(brown)) dataStorage[brown] = std::stod(argMap.at(brown));
            if (!valueLock[enableLinkedCell] && argMap.contains(enableLinkedCell))dataStorage[enableLinkedCell] = std::stoi(argMap.at(enableLinkedCell)) != 0;
            if (!valueLock[enableGrav] && argMap.contains(enableGrav))dataStorage[enableGrav] = std::stoi(argMap.at(enableGrav)) != 0;
            if (!valueLock[enableMembrane] && argMap.contains(enableMembrane))dataStorage[enableMembrane] = std::stoi(argMap.at(enableMembrane)) != 0;
            if (!valueLock[enableMembranePull] && argMap.contains(enableMembranePull))dataStorage[enableMembranePull] = std::stoi(argMap.at(enableMembranePull)) != 0;
            if (!valueLock[enableOMP] && argMap.contains(enableOMP))dataStorage[enableOMP] = std::stoi(argMap.at(enableOMP)) != 0;
            if (!valueLock[rCutoff] && argMap.contains(rCutoff)) dataStorage[rCutoff] = std::stod(argMap.at(rCutoff));
            if (!valueLock[boundingBox_X0] && argMap.contains(boundingBox_X0))dataStorage[boundingBox_X0] = std::stod(argMap.at(boundingBox_X0));
            if (!valueLock[boundingBox_X1] && argMap.contains(boundingBox_X1))dataStorage[boundingBox_X1] = std::stod(argMap.at(boundingBox_X1));
            if (!valueLock[boundingBox_X2] && argMap.contains(boundingBox_X2))dataStorage[boundingBox_X2] = std::stod(argMap.at(boundingBox_X2));
            if (!valueLock[boundCondFront] && argMap.contains(boundCondFront))dataStorage[boundCondFront] = sim::physics::bounds::stot(argMap.at(boundCondFront));
            if (!valueLock[boundCondRear] && argMap.contains(boundCondRear))dataStorage[boundCondRear] = sim::physics::bounds::stot(argMap.at(boundCondRear));
            if (!valueLock[boundCondLeft] && argMap.contains(boundCondLeft))dataStorage[boundCondLeft] = sim::physics::bounds::stot(argMap.at(boundCondLeft));
            if (!valueLock[boundCondRight] && argMap.contains(boundCondRight))dataStorage[boundCondRight] = sim::physics::bounds::stot(argMap.at(boundCondRight));
            if (!valueLock[boundCondTop] && argMap.contains(boundCondTop))dataStorage[boundCondTop] = sim::physics::bounds::stot(argMap.at(boundCondTop));
            if (!valueLock[boundCondBottom] && argMap.contains(boundCondBottom))dataStorage[boundCondBottom] = sim::physics::bounds::stot(argMap.at(boundCondBottom));
            if (!valueLock[dimensions] && argMap.contains(dimensions))dataStorage[dimensions] = std::stoi(argMap.at(dimensions));
            if (!valueLock[logLevel] && argMap.contains(logLevel))dataStorage[logLevel] = std::stoi(argMap.at(logLevel));
            if (!valueLock[benchmark] && argMap.contains(benchmark))dataStorage[benchmark] = std::stoi(argMap.at(benchmark)) != 0;
            if (!valueLock[benchmarkType] && argMap.contains(benchmarkType))dataStorage[benchmarkType] = argMap.at(benchmarkType);
            if (!valueLock[benchMaxBodySize] && argMap.contains(benchMaxBodySize))dataStorage[benchMaxBodySize] = std::stoi(argMap.at(benchMaxBodySize));
            if (!valueLock[benchIterationCount] && argMap.contains(benchIterationCount))dataStorage[benchIterationCount] = std::stoi(argMap.at(benchIterationCount));
            if (!valueLock[enableThermo] && argMap.contains(enableThermo))dataStorage[enableThermo] = std::stoi(argMap.at(enableThermo)) != 0;
            if (!valueLock[thermoTInit] && argMap.contains(thermoTInit))dataStorage[thermoTInit] = std::stod(argMap.at(thermoTInit));
            if (!valueLock[thermoNTerm] && argMap.contains(thermoNTerm))dataStorage[thermoNTerm] = std::stoi(argMap.at(thermoNTerm));
            if (!valueLock[thermoTTarget] && argMap.contains(thermoTTarget))dataStorage[thermoTTarget] = std::stod(argMap.at(thermoTTarget));
            if (!valueLock[thermoDelta_t] && argMap.contains(thermoDelta_t))dataStorage[thermoDelta_t] = std::stod(argMap.at(thermoDelta_t));
            if (!valueLock[enableCheckpointing] && argMap.contains(enableCheckpointing))dataStorage[enableCheckpointing] = std::stoi(argMap.at(enableCheckpointing)) != 0;
            if (!valueLock[gGrav0] && argMap.contains(gGrav0)) dataStorage[gGrav0] = std::stod(argMap.at(gGrav0));
            if (!valueLock[gGrav1] && argMap.contains(gGrav1)) dataStorage[gGrav1] = std::stod(argMap.at(gGrav1));
            if (!valueLock[gGrav2] && argMap.contains(gGrav2)) dataStorage[gGrav2] = std::stod(argMap.at(gGrav2));
            if (!valueLock[simLastIteration] && argMap.contains(simLastIteration))dataStorage[simLastIteration] = std::stoi(argMap.at(simLastIteration));
            if (!valueLock[enableProfiler] && argMap.contains(enableProfiler))dataStorage[enableProfiler] = std::stoi(argMap.at(enableProfiler))!=0;
            if (!valueLock[profilerNumBins] && argMap.contains(profilerNumBins))dataStorage[profilerNumBins] = std::stoi(argMap.at(profilerNumBins));
            if (!valueLock[thermoType_t] && argMap.contains(thermoType_t))dataStorage[thermoType_t] = ((std::stoi(argMap.at(thermoType_t)))==0)?ThermoMode::normalMode:ThermoMode::pipeMode;
        }
    };

} // input
