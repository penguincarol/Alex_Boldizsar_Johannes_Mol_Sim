//
// Created by alex on 30.11.22.
//

#pragma once

#include <unordered_map>
#include <string>
#include "ForceFunctorBase.h"

namespace sim::physics::force {
    enum type {
        gravity,
        lennardJones,
        types_count [[maybe_unused]]
    };

    extern std::unordered_map<std::string, type> type_map;

    /**
     * Converts the given string to a known type.
     * */
    type stot(const std::string &);

} //sim::physics::force
