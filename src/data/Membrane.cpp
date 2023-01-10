//
// Created by johnny on 09.01.23.
//

#include "Membrane.h"

#include<vector>

Membrane::Membrane(double k_r, double r0_r, std::vector<std::tuple<unsigned long, unsigned long>> &springs_r): k(k_r), r0(r0_r), springs(springs_r){
}