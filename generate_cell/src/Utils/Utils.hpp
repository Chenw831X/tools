//
// Created by Wei Chen on 3/2/22
//

#ifndef UTILS_HPP
#define UTILS_HPP

#include <Eigen/Eigen>

namespace SIM{

class Utils{
public:
    static bool writeMatrixXd(const std::string &filePath, const Eigen::MatrixXd &mat);
    static bool writeMatrixXi(const std::string &filePath, const Eigen::MatrixXi &mat);
    static bool readMatrix(const std::string &filePath, int row, int col, Eigen::MatrixXd &mat);

    // write mesh information to '.txt', prepare for meshio
    // V.size(): number of elements in mesh
    // m = 1: V[i] is (1, 3), point mesh
    // m = 2: V[i] is (2, 3), line mesh
    // m = 4: V[i] is (4, 3), quad mesh
    // m = 8: V[i] is (8, 3), hex mesh
    static bool writeMesh(const std::string &filePath, const std::vector<Eigen::MatrixXd> &V, int m);
    static bool writeMesh(const std::string &filePath, const Eigen::MatrixXd &V, const Eigen::MatrixXi &C);
};

} // namespace SIM

#endif // UTILS_HPP