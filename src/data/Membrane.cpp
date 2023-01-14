//
// Created by johnny on 09.01.23.
//

#include "Membrane.h"

#include <utility>
#include<vector>
#include<iostream>

Membrane::Membrane(double k_, double r0_, std::vector<std::vector<unsigned long>> &membrNodes_, double pullEt_, std::array<double, 3> pullF_, std::vector<std::array<size_t, 2>> pullInd)
:k{k_}, r0(r0_), membrNodes(membrNodes_), pullEt(pullEt_), pullF(pullF_), pullIndices(std::move(pullInd)) {}

void Membrane::printMembrNodeStructure() {
    for(int i = 0; i < membrNodes.size(); i++){

        for(int j = 0; j < membrNodes[i].size(); j++){
            std::cout<<membrNodes[i][j];

            if(j+1< membrNodes.size()){std::cout<<"--";}
        }
        std::cout<<std::endl;

        if(i+1 < membrNodes.size()){
            for(int j = 0; j < membrNodes[i].size();j++){
                std::cout<<"|  ";
            }
            std::cout<<std::endl;
        }
    }
}


std::vector<std::vector<unsigned long>>& Membrane::getMembrNodes(){
    return membrNodes;
}

double Membrane::getDesiredDistance() const {
    return r0;
}

double Membrane::getSpringStrength() const {
    return k;
}

double Membrane::getPullEndTime() const  {
    return pullEt;
}

std::array<double, 3> Membrane::getPullForce() const  {
    return pullF;
}

const std::vector<std::array<size_t, 2>> &Membrane::getPullIndices() const  {
    return pullIndices;
}
