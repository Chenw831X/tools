//
// Created by Wei Chen on 3/2/22
//

#include <iostream>
#include <spdlog/spdlog.h>
#include "Config.hpp"

namespace SIM{

Config::Config(){
    // triangle mesh file path of physical domain
    PD_FilePath1 = "../input/triMeshes/bunny.obj";
    // Mesh Setup
    meshDivision1 = Eigen::RowVector3i(10, 10, 10);
}

} // namespace SIM