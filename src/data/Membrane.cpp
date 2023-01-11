//
// Created by johnny on 09.01.23.
//

#include "Membrane.h"

#include<vector>
#include<iostream>

Membrane::Membrane(double k_, double r0_, std::vector<std::vector<unsigned long>> &membrNodes_):k{k_}, r0(r0_), membrNodes(membrNodes_) {}

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


std::vector<std::vector<unsigned long>> Membrane::getMembrNodes(){
    return membrNodes;
};