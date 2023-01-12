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

    double pullEt;
    double pullF;

public:
    Membrane(double k, double r0, std::vector<std::vector<unsigned long>> &membrNodes, double pullEt = 0, double pullF = 0);

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