//
// main.cpp
// objTrans
//
// Created by Wei Chen on 1/15/22
//

#include <iostream>
#include <string>
#include <Eigen/Eigen>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

int main(){
    std::string Input = "/home/cw/MyCode/CGCourseCode1_cpp/input/Puck_on_Toadstool_2048K.obj";
    std::string Output = "/home/cw/MyCode/CGCourseCode1_cpp/input/Puck_on_Toadstool_2048K.obj";
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(Input, V, F);

    V.conservativeResize(V.rows(), 4);
    V.col(3) = Eigen::MatrixXd::Ones(V.rows(), 1);

    glm::mat4 trans(1.0f);
    trans = glm::rotate(trans, glm::radians(180.0f), glm::vec3(1.0, 0.0, 0.0));
    Eigen::Matrix4d _trans;
    for(int i=0; i<4; ++i){
        for(int j=0; j<4; ++j){
            _trans(i, j) = 1.0 * trans[i][j];
        }
    }

    for(int vI=0; vI<V.rows(); ++vI){
        V.row(vI) = (_trans * V.row(vI).transpose()).transpose();
    }

    V.conservativeResize(V.rows(), 3);
    igl::write_triangle_mesh(Output, V, F);
}