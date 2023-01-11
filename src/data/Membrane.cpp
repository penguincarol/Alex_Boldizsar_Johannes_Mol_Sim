//
// Created by johnny on 09.01.23.
//

#include "Membrane.h"

#include<vector>

Membrane::Membrane(double k_, double r0_, std::vector<std::vector<unsigned long>> &&membrNodes_):k{k_}, r0(r0_), membrNodes(membrNodes_) {}