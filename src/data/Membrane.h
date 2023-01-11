//
// Created by johnny on 09.01.23.
//
#pragma once


#include <vector>

class Membrane {
private:
    double k;
    double r0;
    std::vector<std::vector<unsigned long>> membrNodes;

public:
    Membrane(double k, double r0, std::vector<std::vector<unsigned long>>&& membrNodes);
};