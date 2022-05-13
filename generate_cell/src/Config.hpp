//
// Created by Wei Chen on 3/2/22
//

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <Eigen/Eigen>

namespace SIM{

class Config{
public:
    // triangle mesh file path of physical domain
    std::string PD_FilePath1;
    // Mesh Setup
    Eigen::RowVector3i meshDivision1;

public:
    Config();
};

} // namespace SIM

#endif // CONFIG_HPP