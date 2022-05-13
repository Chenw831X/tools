//
// Created by Wei Chen on 3/2/22
//

#ifndef MESH_HPP
#define MESH_HPP

#include <set>
#include <Eigen/Eigen>

#include "PhysicalDomain.hpp"
#include "Element.hpp"

namespace SIM{

class Mesh{
public:
    Eigen::RowVector3d meshOrigin;
    Eigen::RowVector3d meshLength;
    Eigen::RowVector3i meshDivision;
    Eigen::RowVector3d meshStep;

    std::shared_ptr<PhysicalDomain> physicalDomain;
    int nNodeOneDim; // number of nodes in one edge

public:
    int nNode; // number of nodes in this mesh
    int nDof; // number of dofs in this mesh

    // nodes0, edges0, faces0 of original regular mesh
    int nodeNum0, edgeNum0, faceNum0;
    std::vector<Topology*> nodes0;
    std::vector<Topology*> edges0;
    std::vector<Topology*> faces0;

    // nodes, edges, faces of simulation mesh
    // i.e., remove complete void element from original regular mesh
    int nodeNum, edgeNum, faceNum;
    std::vector<Topology*> nodes;
    std::vector<Topology*> edges;
    std::vector<Topology*> faces;

    // mapping from original regular mesh to simulation mesh
    Eigen::VectorXi n02n; // nodes0 id to nodes id, -1: removed
    Eigen::VectorXi e02e; // edges0 id to edges id
    Eigen::VectorXi f02f; // faces0 id to faces id
    Eigen::VectorXi solid2ele; // original solid to element id, -1: void solid

    int nEle; // number of macro elements in this mesh
    int eleNodeNum; // number of nodes in one element
    int eleDofNum; // number of dofs in one element
    std::vector<Element*> elements;

    // location vector for the elements, store nodes id in each element
    std::vector<Eigen::VectorXi> locationVecs;
    // location dofs for the elements, store dofs id in each element
    std::vector<Eigen::VectorXi> locationDofs;

public:
    Mesh(Eigen::RowVector3i p_meshDivision, std::shared_ptr<PhysicalDomain> p_physicalDomain);
    ~Mesh();

public:
    void createNodes0();
    void createEdges0();
    void createFaces0();
    void createElements();
    void computeFeatures();

    // check relation of a given cube and physical domain
    // (cube is represented by 'Coords') 'Coords': (8, 3)
    // 0: empty, 1: partial, 2: full
    int getRelation(const Eigen::Matrix<double, 8, 3> &Coords);
    // for point P with a global coordinate, get its responding element index
    int getEleIDForPoint(const Eigen::RowVector3d &P) const;

};

} // namespace SIM

#endif // MESH_HPP