//
// Created by Wei Chen on 3/4/22
//

#include <spdlog/spdlog.h>

#include "Topology.hpp"

namespace SIM{

Topology::Topology(const Eigen::MatrixXd &p_Coords, int p_nNodeOneDim) :
        Coords(p_Coords), nNodeOneDim(p_nNodeOneDim){

    assert(Coords.cols()==3);
}

int Topology::assignVIndex(int idx){
    vIndex.setZero(nNode);
    dofIndex.setZero(nNode*DIM_);

    int ret = idx+nNode;
    vIndex = Eigen::VectorXi::LinSpaced(nNode, idx, ret-1);
    dofIndex = Eigen::VectorXi::LinSpaced(nNode*DIM_, idx*3, ret*3-1);

    return ret;
}

Node::Node(const Eigen::MatrixXd &p_Coords, int p_nNodeOneDim) :
        Base(p_Coords, p_nNodeOneDim){

    assert(Base::Coords.rows()==1);
    Base::nNode = 1;

    V0.setZero(Base::nNode, DIM_);
    V0.row(0) = Base::Coords.row(0);
}

Edge::Edge(const Eigen::MatrixXd &p_Coords, int p_nNodeOneDim):
        Base(p_Coords, p_nNodeOneDim){

    assert(Base::Coords.rows()==2);
    Base::nNode = Base::nNodeOneDim - 2;

    V0.setZero(Base::nNode, DIM_);
    Eigen::VectorXd LS = Eigen::VectorXd::LinSpaced(Base::nNodeOneDim, 0.0, 1.0);
    Eigen::RowVector3d dir = Base::Coords.row(1) - Base::Coords.row(0);

    for(int vI=0; vI<Base::nNodeOneDim-2; ++vI){
        V0.row(vI) = Base::Coords.row(0) + dir * LS(vI+1);
    }
}

Face::Face(const Eigen::MatrixXd &p_Coords, int p_nNodeOneDim):
        Base(p_Coords, p_nNodeOneDim){

    assert(Base::Coords.rows()==4);
    Base::nNode = (Base::nNodeOneDim-2) * (Base::nNodeOneDim-2);

    V0.setZero(Base::nNode, DIM_);
    Eigen::VectorXd LS = Eigen::VectorXd::LinSpaced(Base::nNodeOneDim, 0.0, 1.0);
    Eigen::RowVector3d dir1 = Base::Coords.row(1) - Base::Coords.row(0);
    Eigen::RowVector3d dir3 = Base::Coords.row(3) - Base::Coords.row(0);

    // first direction 01, then direction 03
    for(int vI=0; vI<Base::nNodeOneDim-2; ++vI){ // direction 03
        for(int vJ=0; vJ<Base::nNodeOneDim-2; ++vJ){ // direction 01
            V0.row(vI*(Base::nNodeOneDim-2)+vJ) = Base::Coords.row(0) + dir3 * LS(vI+1) + dir1 * LS(vJ+1);
        }
    }
}

} // namespace SIM