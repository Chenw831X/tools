//
// Created by Wei Chen on 11/30/21
//

#include <iostream>
#include <utility>
#include <spdlog/spdlog.h>

#include "Type.hpp"
#include "Element.hpp"
#include "Utils.hpp"

namespace SIM{

Element::Element(int p_eleNodeNum, std::vector<Topology*> p_eleNodes, std::vector<Topology*> p_eleEdges,
                 std::vector<Topology*> p_eleFaces): nNode(p_eleNodeNum), eleNodes(std::move(p_eleNodes)),
                 eleEdges(std::move(p_eleEdges)), eleFaces(std::move(p_eleFaces)){

    assert(eleNodes.size()==8);
    assert(eleEdges.size()==12);
    assert(eleFaces.size()==6);
    nDof = nNode * DIM_;
}

void Element::setupLocationVec(Eigen::VectorXi &locationVec, Eigen::VectorXi &locationDof){
    locationVec.resize(nNode);
    locationDof.resize(nDof);
    int cnt = 0;
    
    // Nodes
    for(int nI=0; nI<8; ++nI){
        int num = eleNodes[nI]->nNode;
        locationVec(Eigen::seqN(cnt, num)) = eleNodes[nI]->vIndex;
        locationDof(Eigen::seqN(cnt*3, num*3)) = eleNodes[nI]->dofIndex;
        cnt += num;
    }
    
    // Edges
    for(int eI=0; eI<12; ++eI){
        int num = eleEdges[eI]->nNode;
        locationVec(Eigen::seqN(cnt, num)) = eleEdges[eI]->vIndex;
        locationDof(Eigen::seqN(cnt*3, num*3)) = eleEdges[eI]->dofIndex;
        cnt += num;
    }

    // Faces
    for(int fI=0; fI<6; ++fI){
        int num = eleFaces[fI]->nNode;
        locationVec(Eigen::seqN(cnt, num)) = eleFaces[fI]->vIndex;
        locationDof(Eigen::seqN(cnt*3, num*3)) = eleFaces[fI]->dofIndex;
        cnt += num;
    }

    // std::cout << "cnt: " << cnt << std::endl;
    assert(cnt == nNode);
}

} // namespace SIM