//
// Created by alex on 30.11.22.
//

#include "types.h"
#include "io/input/cli/CLIArgsParser.h"
#include "io/output/Logging.h"

namespace sim::physics::force {
    std::unordered_map<std::string, type> type_map = {{"gravity",    gravity},
                                                      {"lennardjones", lennardJones}};

    type stot(const std::string &str) {
        auto lowercase = [](std::string &str) {
            std::transform(str.begin(), str.end(), str.begin(), [](unsigned char c) { return std::tolower(c); });
        };
        std::string key = str;
        lowercase(key);
        if (!type_map.contains(key)) io::input::exitFormatError("force type unknown!");
        else return type_map[key];
    }

} // sim::physics::force

