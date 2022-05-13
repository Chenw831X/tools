//
// Created by Wei Chen on 3/1/22
//

#include <iostream>
#include <chrono>
#include <Eigen/Eigen>
#include <spdlog/spdlog.h>

#include "Type.hpp"
#include "Config.hpp"
#include "Mesh.hpp"

SIM::Config config;
std::shared_ptr<SIM::PhysicalDomain> physicalDomain;
std::shared_ptr<SIM::Mesh> mesh;

int main(int argc, char **argv){
    auto programStart = std::chrono::steady_clock::now();

    // 0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical, 6=off
    spdlog::set_level(spdlog::level::trace); // print all logs
    Eigen::setNbThreads(8);
    spdlog::info("OMP/Eigen: number of threads: {}", Eigen::nbThreads());


    physicalDomain = std::make_shared<SIM::PhysicalDomain>(config.PD_FilePath1);
    mesh = std::make_shared<SIM::Mesh>(config.meshDivision1, physicalDomain);


    std::chrono::duration<double> program_(std::chrono::steady_clock::now() - programStart);
    spdlog::info("program end, total time cost: {:.6f} s", program_.count());
    return 0;
}