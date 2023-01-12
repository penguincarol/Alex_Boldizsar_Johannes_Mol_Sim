//
// Created by johnny on 09.01.23.
//
#pragma once


#include <vector>
#include <array>

class Membrane {
private:
    double k;
    double r0;
    std::vector<std::vector<unsigned long>> membrNodes;

    double pullEt;
    std::array<double,3> pullF;
    std::vector<std::array<size_t,2>> pullIndices;

public:
    Membrane(double k, double r0, std::vector<std::vector<unsigned long>> &membrNodes, double pullEt = 0, std::array<double,3> pullF = {0}, std::vector<std::array<size_t, 2>> pullInd = {});

    /**
     * Helper method that prints out the current membrane node structure
     * For debugging purposes only
     * only really makes sense for small membranes
     */
    void printMembrNodeStructure();

    /**
     * not supposed to be used in performance critical areas due to unnessecary copying
     * @return
     */
    std::vector<std::vector<unsigned long>>& getMembrNodes();

    double getDesiredDistance();

    double getSpringStrength();
};