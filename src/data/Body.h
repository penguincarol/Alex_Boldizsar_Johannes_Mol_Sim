#pragma once

#include "Eigen"

enum Shape {
    cuboid, sphere, particle, membrane, num_shapes      //Alex you sick little bastard, that is the sickest abuse of enums as int i have ever seen
};
const std::array<enum Shape, num_shapes> all_shapes = {cuboid, sphere, particle, membrane};

struct Body {
    Shape shape;    /**<defines the shape of the body and therefore the interpretation of the following values*/
    Eigen::Vector3d fixpoint; /**<front bottom left point of cuboid, middle of sphere, etc*/
    Eigen::Vector3d dimensions; /**<amount of particles in each dimensional direction; Sphere assumes that dimensions[0] == dimensions[1] && dimensions[1] == dimensions[2]*/
    double distance;
    double mass;
    Eigen::Vector3d start_velocity;

    double desiredDistance; /**These two parameters only get used in case that shape is Membrane. Otherwise it's undefined*/
    double springStrength; /**These two parameters only get used in case that shape is Membrane. Otherwise it's undefined*/
    double pullEndTime;
    Eigen::Vector3d pullForce;
    std::vector<std::array<size_t, 2>> pullIndices;
};