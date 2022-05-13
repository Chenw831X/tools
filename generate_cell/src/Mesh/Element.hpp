//
// Created by Wei Chen on 11/29/21
//

#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <Eigen/Eigen>

#include "Topology.hpp"

namespace SIM{

class Element{
public:
    std::vector<Topology*> eleNodes;
    std::vector<Topology*> eleEdges;
    std::vector<Topology*> eleFaces;

    int nNode; // total node (not Node object) in element
    int nDof; // total dof in element (nNode * DIM)

public:
    Element(int p_eleNodeNum, std::vector<Topology*> p_eleNodes,
            std::vector<Topology*> p_eleEdges, std::vector<Topology*> p_eleFaces);

public:
    void setupLocationVec(Eigen::VectorXi &locationVec, Eigen::VectorXi &locationDof);

};

} // namespace

#endif // ELEMENT_HPP