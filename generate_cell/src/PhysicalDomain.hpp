//
// Created by Wei Chen on 3/2/22
//

#ifndef PHYSICALDOMAIN_HPP
#define PHYSICALDOMAIN_HPP

#include <Eigen/Eigen>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/fast_winding_number.h>

namespace SIM{

class PhysicalDomain{
private:
    std::string filePath;
    igl::FastWindingNumberBVH fwn_bvh;

public:
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    int numV, numF;

    Eigen::Matrix<double, 2, 3> bbox; // bounding box of physical domain
    Eigen::RowVector3d lenBBox; // length of bounding box in 3 dimensions

public:
    explicit PhysicalDomain(std::string p_filePath);

public:
    // given query Q (global coordinates), return W (their domain index)
    // 0: fictitious domain
    // 1: physical domain
    // domain index corresponds to the material index
    void getDomainID(const Eigen::MatrixXd &Q, Eigen::VectorXi &domainID) const;
};

} // namespace SIM

#endif // PHYSICALDOMAIN_HPP