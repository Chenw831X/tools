//
// Created by Wei Chen on 3/2/22
//

#include <iostream>
#include <utility>
#include <spdlog/spdlog.h>

#include "PhysicalDomain.hpp"

namespace SIM{

PhysicalDomain::PhysicalDomain(std::string p_filePath): filePath(std::move(p_filePath)){
    if(!igl::read_triangle_mesh(filePath, V, F)){
        spdlog::error("error on reading triangle mesh from " + filePath);
        exit(-1);
    }
    assert(V.cols() == 3);
    assert(F.cols() == 3);

    numV = (int)V.rows();
    numF = (int)F.rows();

    bbox.row(0) = V.colwise().minCoeff();
    bbox.row(1) = V.colwise().maxCoeff();
    lenBBox = bbox.row(1) - bbox.row(0);

    igl::fast_winding_number(V, F, 2, fwn_bvh);

    spdlog::info("physical domain constructed(" + filePath + "): " +
        std::to_string(numV) + " points, " + std::to_string(numF) + " triangles");
}

void PhysicalDomain::getDomainID(const Eigen::MatrixXd &Q, Eigen::VectorXi &domainID) const{
    assert(Q.cols()==3);
    // winding number bigger than val will be seen inside physical domain, otherwise outside
    double val = 0.5;

    Eigen::VectorXd W;
    igl::fast_winding_number(fwn_bvh, 2, Q, W);
    assert(W.rows() == Q.rows());

    domainID = (W.array() > val).cast<int>();
}
    
} // namespace SIM