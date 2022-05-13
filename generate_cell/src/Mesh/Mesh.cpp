//
// Created by Wei Chen on 3/2/22
//

#include <iostream>
#include <utility>
#include <spdlog/spdlog.h>
#include <igl/write_triangle_mesh.h>

#include "Type.hpp"
#include "Utils.hpp"
#include "Mesh.hpp"

namespace SIM{

Mesh::Mesh(Eigen::RowVector3i p_meshDivision, std::shared_ptr<PhysicalDomain> p_physicalDomain):
        meshDivision(std::move(p_meshDivision)), physicalDomain(std::move(p_physicalDomain)){

    meshOrigin = physicalDomain->bbox.row(0) - physicalDomain->lenBBox * 0.05;
    meshLength = physicalDomain->lenBBox * 1.10;
    meshStep = meshLength.array() / meshDivision.cast<double>().array();

    nNodeOneDim = 2;
    eleNodeNum = 8;
    eleDofNum = eleNodeNum * DIM_;

    createNodes0();
    createEdges0();
    createFaces0();
    createElements();

    nNode = 0;
    Eigen::MatrixXd mesh_V(0, 3);
    for(int nI=0; nI<nodeNum; ++nI){
        int tmp = nNode;
        nNode = nodes[nI]->assignVIndex(nNode);
        mesh_V.conservativeResize(nNode, 3);
        mesh_V(Eigen::seq(tmp, Eigen::last), Eigen::all) = nodes[nI]->Coords;
    }
    for(int eI=0; eI<edgeNum; ++eI){
        int tmp = nNode;
        nNode = edges[eI]->assignVIndex(nNode);
        mesh_V.conservativeResize(nNode, 3);
        mesh_V(Eigen::seq(tmp, Eigen::last), Eigen::all) = edges[eI]->Coords;
    }
    for(int fI=0; fI<faceNum; ++fI){
        int tmp = nNode;
        nNode = faces[fI]->assignVIndex(nNode);
        mesh_V.conservativeResize(nNode, 3);
        mesh_V(Eigen::seq(tmp, Eigen::last), Eigen::all) = faces[fI]->Coords;
    }
    nDof = nNode * DIM_;
    computeFeatures();

    Eigen::MatrixXi mesh_C(nEle, eleNodeNum);
    for(int eleI=0; eleI<nEle; ++eleI){
        mesh_C.row(eleI) = locationVecs[eleI].transpose();
    }

    Utils::writeMesh("/home/cw/MyCode/tools/generate_cell/output/mesh.txt", mesh_V, mesh_C);
    Utils::writeMatrixXd("/home/cw/MyCode/tools/generate_cell/output/mesh_V.txt", mesh_V);
    Utils::writeMatrixXi("/home/cw/MyCode/tools/generate_cell/output/mesh_C.txt", mesh_C);

    spdlog::info("original regular mesh constructed: from ({}, {}, {}) to ({}, {}, {}), with ({}, {}, {}) solids",
        meshOrigin(0), meshOrigin(1), meshOrigin(2), meshOrigin(0)+meshLength(0),
        meshOrigin(1)+meshLength(1), meshOrigin(2)+meshLength(2), meshDivision(0),
        meshDivision(1), meshDivision(2));
    spdlog::info("number of nodes: {}, number of dofs: {} in simulation mesh", nNode, nDof);
    spdlog::info("number of elements in simulation mesh: {}", nEle);

}

Mesh::~Mesh(){
    for(int nI=0; nI < nodeNum0; ++nI){
        delete nodes0[nI];
    }

    for(int eI=0; eI < edgeNum0; ++eI){
        delete edges0[eI];
    }

    for(int fI=0; fI < faceNum0; ++fI){
        delete faces0[fI];
    }

    for(int eleI=0; eleI<nEle; ++eleI){
        delete elements[eleI];
    }

}

void Mesh::createNodes0(){
    spdlog::info("create nodes0 ...");

    bool debug = false;
    // bool debug = true;

    nodeNum0 = (meshDivision(0) + 1) * (meshDivision(1) + 1) * (meshDivision(2) + 1);
    nodes0.resize(nodeNum0);

    double xIncrement = meshLength(0) / (1.0 * meshDivision(0));
    double yIncrement = meshLength(1) / (1.0 * meshDivision(1));
    double zIncrement = meshLength(2) / (1.0 * meshDivision(2));
    
    double xx=meshOrigin(0), yy=meshOrigin(1), zz=meshOrigin(2);
    int cnt = 0;
    for(int k=0; k<=meshDivision(2); ++k){
        for(int j=0; j<=meshDivision(1); ++j){
            for(int i=0; i<=meshDivision(0); ++i){
                if(debug){
                    std::cout << "Node #" << cnt << ":  " << "(" << xx <<
                        ", " << yy << ", " << zz << ")" << std::endl;
                }
                nodes0[cnt] = new Node(Eigen::Matrix<double, 1, 3>(xx, yy, zz), nNodeOneDim);
                ++cnt;

                xx += xIncrement;
            }

            xx = meshOrigin(0);
            yy += yIncrement;
        }

        yy = meshOrigin(1);
        zz += zIncrement;
    }
    assert(cnt == nodeNum0);

    // print for debug
    if(debug){
        std::cout << "cnt: " << cnt << std::endl;
        std::cout << "nodeNum0: " << nodeNum0 << std::endl;
    }
}

void Mesh::createEdges0(){
    spdlog::info("create edges0 ...");

    bool debug = false;
    // bool debug = true;

    edgeNum0 = (meshDivision(2) + 1) * (meshDivision(1) + 1) * meshDivision(0) +
        (meshDivision(2) + 1) * (meshDivision(0) + 1) * meshDivision(1) +
        (meshDivision(1) + 1) * (meshDivision(0) + 1) * meshDivision(2);
    edges0.resize(edgeNum0);

    int cnt = 0;

    // create edges in x direction
    int nNodesInXY = (meshDivision(0) + 1) * (meshDivision(1) + 1);
    std::vector<int> nodesInXY(nNodesInXY);
    for(int i=0; i<nNodesInXY; ++i){
        nodesInXY[i] = i;
    }

    // outest loop in z direction
    for(int k=0; k<=meshDivision(2); ++k){
        for(int j=0; j<=meshDivision(1); ++j){
            int startI = j * (meshDivision(0) + 1);
            for(int i=0; i<meshDivision(0); ++i){
                if(debug){
                    std::cout << "Edge #" << cnt << ":  " << nodesInXY[startI+i] <<
                        " -- " << nodesInXY[startI+i+1] << std::endl;
                }
                Eigen::Matrix<double, 2, 3> Coords;
                Coords.row(0) = nodes0[nodesInXY[startI + i]]->Coords.row(0);
                Coords.row(1) = nodes0[nodesInXY[startI + i + 1]]->Coords.row(0);
                edges0[cnt] = new Edge(Coords, nNodeOneDim);
                ++cnt;
            }
        }

        for(int t=0; t<nNodesInXY; ++t){
            nodesInXY[t] += nNodesInXY;
        }
    }


    // create edges in y direction
    int nNodesInYZ = (meshDivision(1) + 1) * (meshDivision(2) + 1);
    std::vector<int> nodesInYZ(nNodesInYZ);
    nodesInYZ[0] = 0;
    for(int i=1; i<=meshDivision(1); ++i){
        nodesInYZ[i] = nodesInYZ[i-1] + meshDivision(0) + 1;
    }
    for(int i=1; i<=meshDivision(2); ++i){
        for(int j=0; j<=meshDivision(1); ++j){
            nodesInYZ[i*(meshDivision(1)+1)+j] = nodesInYZ[(i-1)*(meshDivision(1)+1)+j] +
                nNodesInXY;
        }
    }

    // outest loop in x direction
    for(int i=0; i<=meshDivision(0); ++i){
        for(int k=0; k<=meshDivision(2); ++k){
            int startI = k * (meshDivision(1) + 1);
            for(int j=0; j<meshDivision(1); ++j){
                if(debug){
                    std::cout << "Edge #" << cnt << ":  " << nodesInYZ[startI+j] <<
                        " -- " << nodesInYZ[startI+j+1] << std::endl;
                }
                Eigen::Matrix<double, 2, 3> Coords;
                Coords.row(0) = nodes0[nodesInYZ[startI + j]]->Coords.row(0);
                Coords.row(1) = nodes0[nodesInYZ[startI + j + 1]]->Coords.row(0);
                edges0[cnt] = new Edge(Coords, nNodeOneDim);
                ++cnt;
            }
        }

        for(int t=0; t<nNodesInYZ; ++t){
            nodesInYZ[t] += 1;
        }
    }


    // create edges in z direction
    int nNodesInZX = (meshDivision(0) + 1) * (meshDivision(2) + 1);
    std::vector<int> nodesInZX(nNodesInZX);
    nodesInZX[0] = 0;
    for(int i=1; i<=meshDivision(2); ++i){
        nodesInZX[i] = nodesInZX[i-1] + nNodesInXY;
    }
    for(int i=1; i<=meshDivision(0); ++i){
        for(int j=0; j<=meshDivision(2); ++j){
            nodesInZX[i*(meshDivision(2)+1)+j] = nodesInZX[(i-1)*(meshDivision(2)+1)+j] + 1;
        }
    }

    // outest loop in y direction
    for(int j=0; j<=meshDivision(1); ++j){
        for(int i=0; i<=meshDivision(0); ++i){
            int startI = i * (meshDivision(2) + 1);
            for(int k=0; k<meshDivision(2); ++k){
                if(debug){
                    std::cout << "Edge #" << cnt << ":  " << nodesInZX[startI+k] <<
                        " -- " << nodesInZX[startI+k+1] << std::endl;
                }
                Eigen::Matrix<double, 2, 3> Coords;
                Coords.row(0) = nodes0[nodesInZX[startI + k]]->Coords.row(0);
                Coords.row(1) = nodes0[nodesInZX[startI + k + 1]]->Coords.row(0);
                edges0[cnt] = new Edge(Coords, nNodeOneDim);
                ++cnt;
            }
        }

        for(int t=0; t<nNodesInZX; ++t){
            nodesInZX[t] += meshDivision[0] + 1;
        }
    }

    assert(cnt == edgeNum0);

    // print for debug
    if(debug){
        std::cout << "cnt: " << cnt << std::endl;
        std::cout << "edgeNum: " << edgeNum0 << std::endl;
    }
}

void Mesh::createFaces0(){
    spdlog::info("create faces0 ...");

    bool debug = false;
    // bool debug = true;

    faceNum0 = (meshDivision(2) + 1) * meshDivision(1) * meshDivision(0) +
        (meshDivision(0) + 1) * meshDivision(2) * meshDivision(1) +
        (meshDivision(1) + 1) * meshDivision(2) * meshDivision(0);
    faces0.resize(faceNum0);

    int cnt = 0;

    // create faces in xy area
    int nNodesInXY = (meshDivision(0) + 1) * (meshDivision(1) + 1);
    std::vector<int> nodesInXY(nNodesInXY);
    for(int i=0; i<nNodesInXY; ++i){
        nodesInXY[i] = i;
    }

    // outest loop in z direction
    for(int k=0; k<=meshDivision(2); ++k){
        for(int j=0; j<meshDivision(1); ++j){
            int startI1 = j * (meshDivision(0) + 1);
            int startI2 = (j + 1) * (meshDivision(0) + 1);
            for(int i=0; i<meshDivision(0); ++i){
                if(debug){
                    std::cout << "Face #" << cnt << ": " << nodesInXY[startI1+i] << " -- " <<
                        nodesInXY[startI1+i+1] << " -- " << nodesInXY[startI2+i+1] <<
                        " -- " << nodesInXY[startI2+i] << std::endl;
                }
                Eigen::Matrix<double, 4, 3> Coords;
                Coords.row(0) = nodes0[nodesInXY[startI1 + i]]->Coords.row(0);
                Coords.row(1) = nodes0[nodesInXY[startI1 + i + 1]]->Coords.row(0);
                Coords.row(2) = nodes0[nodesInXY[startI2 + i + 1]]->Coords.row(0);
                Coords.row(3) = nodes0[nodesInXY[startI2 + i]]->Coords.row(0);
                faces0[cnt] = new Face(Coords, nNodeOneDim);
                ++cnt;
            }
        }

        for(int t=0; t<nNodesInXY; ++t){
            nodesInXY[t] += nNodesInXY;
        }
    }


    // create faces in yz area
    int nNodesInYZ = (meshDivision(1) + 1) * (meshDivision(2) + 1);
    std::vector<int> nodesInYZ(nNodesInYZ);
    nodesInYZ[0] = 0;
    for(int i=1; i<=meshDivision(1); ++i){
        nodesInYZ[i] = nodesInYZ[i-1] + meshDivision(0) + 1;
    }
    for(int i=1; i<=meshDivision(2); ++i){
        for(int j=0; j<=meshDivision(1); ++j){
            nodesInYZ[i*(meshDivision(1)+1)+j] = nodesInYZ[(i-1)*(meshDivision(1)+1)+j] +
                nNodesInXY;
        }
    }

    // outest loop in x direction
    for(int i=0; i<=meshDivision(0); ++i){
        for(int k=0; k<meshDivision(2); ++k){
            int startI1 = k * (meshDivision[1] + 1);
            int startI2 = (k + 1) * (meshDivision[1] + 1);
            for(int j=0; j<meshDivision(1); ++j){
                if(debug){
                    std::cout << "Face #" << cnt << ": " << nodesInYZ[startI1+j] << " -- " <<
                        nodesInYZ[startI1+j+1] << " -- " << nodesInYZ[startI2+j+1] <<
                        " -- " << nodesInYZ[startI2+j] << std::endl;
                }
                Eigen::Matrix<double, 4, 3> Coords;
                Coords.row(0) = nodes0[nodesInYZ[startI1 + j]]->Coords.row(0);
                Coords.row(1) = nodes0[nodesInYZ[startI1 + j + 1]]->Coords.row(0);
                Coords.row(2) = nodes0[nodesInYZ[startI2 + j + 1]]->Coords.row(0);
                Coords.row(3) = nodes0[nodesInYZ[startI2 + j]]->Coords.row(0);
                faces0[cnt] = new Face(Coords, nNodeOneDim);
                ++cnt;
            }
        }

        for(int t=0; t<nNodesInYZ; ++t){
            nodesInYZ[t] += 1;
        }
    }


    // create faces in zx area
    int nNodesInZX = (meshDivision(0) + 1) * (meshDivision(2) + 1);
    std::vector<int> nodesInZX(nNodesInZX);
    nodesInZX[0] = 0;
    for(int i=1; i<=meshDivision(2); ++i){
        nodesInZX[i] = nodesInZX[i-1] + nNodesInXY;
    }
    for(int i=1; i<=meshDivision(0); ++i){
        for(int j=0; j<=meshDivision(2); ++j){
            nodesInZX[i*(meshDivision(2)+1)+j] = nodesInZX[(i-1)*(meshDivision(2)+1)+j] + 1;
        }
    }

    // outest loop in direction y
    for(int j=0; j<=meshDivision(1); ++j){
        for(int i=0; i<meshDivision(0); ++i){
            int startI1 = i * (meshDivision(2) + 1);
            int startI2 = (i + 1) * (meshDivision(2) + 1);
            for(int k=0; k<meshDivision(2); ++k){
                if(debug){
                    std::cout << "Face #" << cnt << ": " << nodesInZX[startI1+k] << " -- " <<
                        nodesInZX[startI1+k+1] << " -- " << nodesInZX[startI2+k+1] <<
                        " -- " << nodesInZX[startI2+k] << std::endl;
                }
                Eigen::Matrix<double, 4, 3> Coords;
                Coords.row(0) = nodes0[nodesInZX[startI1 + k]]->Coords.row(0);
                Coords.row(1) = nodes0[nodesInZX[startI1 + k + 1]]->Coords.row(0);
                Coords.row(2) = nodes0[nodesInZX[startI2 + k + 1]]->Coords.row(0);
                Coords.row(3) = nodes0[nodesInZX[startI2 + k]]->Coords.row(0);
                faces0[cnt] = new Face(Coords, nNodeOneDim);
                ++cnt;
            }
        }

        for(int t=0; t<nNodesInZX; ++t){
            nodesInZX[t] += meshDivision(0) + 1;
        }
    }

    assert(cnt == faceNum0);

    // print for debug
    if(debug){
        std::cout << "cnt: " << cnt << std::endl;
        std::cout << "faceNum: " << faceNum0 << std::endl;
    }
}

void Mesh::createElements(){
    spdlog::info("create elements ...");

    bool debug = false;
    bool writeEle = false; // write element mesh to 'mesh_ele.txt'
    std::vector<Eigen::MatrixXd> mesh_ele;

    n02n = Eigen::VectorXi::Constant(nodeNum0, -1);
    e02e = Eigen::VectorXi::Constant(edgeNum0, -1);
    f02f = Eigen::VectorXi::Constant(faceNum0, -1);
    int solidNum0 = meshDivision(2) * meshDivision(1) * meshDivision(0);
    solid2ele = Eigen::VectorXi::Constant(solidNum0, -1);

    int ncnt = 0;
    int ecnt = 0;
    int fcnt = 0;
    int cnt = 0; // element count

    std::vector<int> nodeID(8);
    std::vector<int> edgeID(12);
    std::vector<int> faceID(6);
    std::vector<Topology*> eleNodes(8);
    std::vector<Topology*> eleEdges(12);
    std::vector<Topology*> eleFaces(6);

    int nNodesInXY = (meshDivision(0) + 1) * (meshDivision(1) + 1);
    int nEdgesInX = meshDivision(0) * (meshDivision(1) + 1) * (meshDivision(2) + 1);
    int nEdgesInY = meshDivision(1) * (meshDivision(0) + 1) * (meshDivision(2) + 1);
    int nFacesInXY = (meshDivision(2) + 1) * meshDivision(1) * meshDivision(0); 
    int nFacesInYZ = (meshDivision(0) + 1) * meshDivision(2) * meshDivision(1);

    for(int k=0; k<meshDivision(2); ++k){
        int nStartK = k * nNodesInXY;
        for(int j=0; j<meshDivision(1); ++j){
            int nStartJ = nStartK + j * (meshDivision(0) + 1);
            for(int i=0; i<meshDivision(0); ++i){

                nodeID[0] = nStartJ + i;
                nodeID[1] = nodeID[0] + 1;
                nodeID[3] = nodeID[0] + meshDivision(0) + 1;
                nodeID[2] = nodeID[3] + 1;
                nodeID[4] = nodeID[0] + nNodesInXY;
                nodeID[5] = nodeID[1] + nNodesInXY;
                nodeID[6] = nodeID[2] + nNodesInXY;
                nodeID[7] = nodeID[3] + nNodesInXY;

                // get relation of current solid with model
                Eigen::Matrix<double, 8, 3> solid;
                for(int nI=0; nI<8; ++nI){
                    solid.row(nI) = nodes0[nodeID[nI]]->Coords;
                }
                // remove empty cell
                if(getRelation(solid) == 0){
                    continue;
                }

                if(writeEle){
                    mesh_ele.emplace_back(solid);
                }

                edgeID[0] = k * meshDivision(0) * (meshDivision(1) + 1) +
                    j * meshDivision(0) + i;
                edgeID[2] = edgeID[0] + meshDivision(0);
                edgeID[3] = nEdgesInX + i * meshDivision(1) * (meshDivision(2) + 1) +
                    k * meshDivision(1) + j;
                edgeID[1] = edgeID[3] + meshDivision(1) * (meshDivision(2) + 1);
                edgeID[4] = nEdgesInX + nEdgesInY + j * meshDivision(2) *
                    (meshDivision(0) + 1) + i * meshDivision(2) + k;
                edgeID[5] = edgeID[4] + meshDivision(2);
                edgeID[7] = edgeID[4] + meshDivision(2) * (meshDivision(0) + 1);
                edgeID[6] = edgeID[7] + meshDivision(2);
                edgeID[8] = edgeID[0] + meshDivision(0) * (meshDivision(1) + 1);
                edgeID[9] = edgeID[1] + meshDivision(1);
                edgeID[10] = edgeID[2] + meshDivision(0) * (meshDivision(1) + 1);
                edgeID[11] = edgeID[3] + meshDivision(1);

                faceID[0] = k * meshDivision(0) * meshDivision(1) + j * meshDivision(0) + i;
                faceID[5] = faceID[0] + meshDivision(0) * meshDivision(1);
                faceID[4] = nFacesInXY + i * meshDivision(1) * meshDivision(2) +
                    k * meshDivision(1) + j;
                faceID[2] = faceID[4] + meshDivision(1) * meshDivision(2);
                faceID[1] = nFacesInXY + nFacesInYZ + j * meshDivision(2) * meshDivision(0) +
                    i * meshDivision(2) + k;
                faceID[3] = faceID[1] + meshDivision(2) * meshDivision(0);

                for(int nI=0; nI<8; ++nI){
                    int id0 = nodeID[nI];
                    if(n02n(id0) == -1){
                        n02n(id0) = ncnt;
                        nodes.emplace_back(nodes0[id0]);
                        ++ncnt;
                    }
                    eleNodes[nI] = nodes[n02n(id0)];
                }
                for(int eI=0; eI<12; ++eI){
                    int id0 = edgeID[eI];
                    if(e02e(id0) == -1){
                        e02e(id0) = ecnt;
                        edges.emplace_back(edges0[id0]);
                        ++ecnt;
                    }
                    eleEdges[eI] = edges[e02e(id0)];
                }
                for(int fI=0; fI<6; ++fI){
                    int id0 = faceID[fI];
                    if(f02f(id0) == -1){
                        f02f(id0) = fcnt;
                        faces.emplace_back(faces0[id0]);
                        ++fcnt;
                    }
                    eleFaces[fI] = faces[f02f(id0)];
                }

                int solidID = k * meshDivision(1) * meshDivision(0) + j * meshDivision(0) + i;
                solid2ele(solidID) = cnt;

                if(debug){
                    std::cout << "***** Element #" << cnt << " *****" << std::endl;
                    
                    std::cout << "Node:  ";
                    for(int nI=0; nI<8; ++nI){
                        std::cout << nodeID[nI] << ", ";
                    }
                    std::cout << std::endl;

                    std::cout << "Edge:  ";
                    for(int eI=0; eI<12; ++eI){
                        std::cout << edgeID[eI] << ", ";
                    }
                    std::cout << std::endl;

                    std::cout << "Face:  ";
                    for(int fI=0; fI<6; ++fI){
                        std::cout << faceID[fI] << ", ";
                    }
                    std::cout << std::endl;

                    std::cout << std::endl;
                }
                elements.emplace_back(new Element(eleNodeNum, eleNodes,
                                                  eleEdges, eleFaces));
                ++cnt;
            }
        }
    }

    // setup number of nodes, edges, faces, solids, elements of simulation mesh
    nodeNum = ncnt;
    edgeNum = ecnt;
    faceNum = fcnt;
    nEle = cnt;

    // for(int eleI=0; eleI<nEle; ++eleI){
    //     std::cout << "***** Element #" << eleI << " *****" << std::endl;
    //     elements[eleI]->debug();
    //     std::cout << std::endl;
    // }

    if(writeEle){
        spdlog::info("write element mesh to \'mesh_ele.txt\'");
        Utils::writeMesh("../output/mesh_ele.txt", mesh_ele, 8);
    }

    // print for debug
    if(debug){
        std::cout << "cnt: " << cnt << std::endl;
        std::cout << "nEle: " << nEle << std::endl;
    }
}

// compute locationVecs
void Mesh::computeFeatures(){
    // setup locationVecs
    locationVecs.resize(nEle);
    locationDofs.resize(nEle);
    for(int eleI=0; eleI<nEle; ++eleI){
        elements[eleI]->setupLocationVec(locationVecs[eleI], locationDofs[eleI]);
    }
}

int Mesh::getRelation(const Eigen::Matrix<double, 8, 3> &Coords){
    // sample point method
    // int m = 24; // number of seed points in one dimension
    int m = 8; // number of seed points in one dimension
    Eigen::VectorXd Xseeds = Eigen::VectorXd::LinSpaced(m, Coords(0, 0), Coords(1, 0));
    Eigen::VectorXd Yseeds = Eigen::VectorXd::LinSpaced(m, Coords(0, 1), Coords(3, 1));
    Eigen::VectorXd Zseeds = Eigen::VectorXd::LinSpaced(m, Coords(0, 2), Coords(4, 2));

    int numPoint = m * m * m;
    Eigen::MatrixXd seedPoints(numPoint, 3);
    Eigen::VectorXi W;
    int cnt = 0;
    for(int i=0; i<m; ++i){
        for(int j=0; j<m; ++j){
            for(int k=0; k<m; ++k){
                seedPoints.row(cnt) = Eigen::RowVector3d(Xseeds(i), Yseeds(j), Zseeds(k));
                ++cnt;
            }
        }
    }

    physicalDomain->getDomainID(seedPoints, W);
    int num = W.count();
    // std::cout << "num: " << num << std::endl;
    if(num == 0){
        return 0; // empty
    }
    else if(num == numPoint){
        return 2; // full
    }
    else{
        return 1; // partial
    }
}

int Mesh::getEleIDForPoint(const Eigen::RowVector3d &P) const{
    Eigen::RowVector3i tmp = ((P - meshOrigin).array() / meshStep.array()).cast<int>();
    if(tmp(0) == meshDivision(0)){
        --tmp(0);
    }
    if(tmp(1) == meshDivision(1)){
        --tmp(1);
    }
    if(tmp(2) == meshDivision(2)){
        --tmp(2);
    }
    assert(tmp(0) >= 0 && tmp(0) < meshDivision(0) &&
           tmp(1) >= 0 && tmp(1) < meshDivision(1) &&
           tmp(2) >= 0 && tmp(2) < meshDivision(2));

    int id0 = tmp(2) * meshDivision(0) * meshDivision(1) + tmp(1) * meshDivision(0) + tmp(0);
    int ret = solid2ele(id0);
    assert(ret >=0 && ret < nEle);
    return ret;
}

} // namespace SIM