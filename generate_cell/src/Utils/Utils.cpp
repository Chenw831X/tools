//
// Created by Wei Chen on 3/2/22
//

#include "iostream"
#include <cmath>
#include <spdlog/spdlog.h>

#include <time.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include "Utils.hpp"

namespace SIM{

bool Utils::writeMatrixXd(const std::string &filePath, const Eigen::MatrixXd &mat){
    FILE* out = fopen(filePath.c_str(), "w");
    if(!out){
        spdlog::error("error on writing matrix");
        return false;
    }

    int row = (int)mat.rows(), col = (int)mat.cols();
    for(int i=0; i<row; ++i){
        for(int j=0; j<col; ++j){
            if(j){
                fprintf(out, " ");
            }
            fprintf(out, "%.18le", mat(i, j));
        }
        fprintf(out, "\n");
    }

    fclose(out);
    return true;
}

bool Utils::writeMatrixXi(const std::string &filePath, const Eigen::MatrixXi &mat){
    FILE* out = fopen(filePath.c_str(), "w");
    if(!out){
        spdlog::error("error on writing matrix");
        return false;
    }

    int row = (int)mat.rows(), col = (int)mat.cols();
    for(int i=0; i<row; ++i){
        for(int j=0; j<col; ++j){
            if(j){
                fprintf(out, " ");
            }
            fprintf(out, "%d", mat(i, j));
        }
        fprintf(out, "\n");
    }

    fclose(out);
    return true;
}

bool Utils::readMatrix(const std::string &filePath, int row, int col, Eigen::MatrixXd &mat){
    FILE* in = fopen(filePath.c_str(), "r");
    if(!in){
        spdlog::error("error on reading matrix");
        return false;
    }

    mat.resize(row, col);
    for(int i=0; i<row; ++i){
        for(int j=0; j<col; ++j){
            fscanf(in, "%le", &mat(i, j));
        }
    }

    fclose(in);
    return true;
}

bool Utils::writeMesh(const std::string &filePath, const std::vector<Eigen::MatrixXd> &V,
                      int m){
    FILE* out = fopen(filePath.c_str(), "w");
    if(!out){
        std::cout << "error on writing mesh" << std::endl;
        return false;
    }

    int nEle = V.size();
    // write number of V and F
    fprintf(out, "%d\n%d\n", m*nEle, nEle);

    for(int eleI=0; eleI<nEle; ++eleI){
        assert(m == V[eleI].rows());
        for(int i=0; i<m; ++i){
            fprintf(out, "%le %le %le\n", V[eleI](i, 0), V[eleI](i, 1), V[eleI](i, 2));
        }
    }
    int cnt = 0;
    for(int eleI=0; eleI<nEle; ++eleI){
        for(int i=0; i<m; ++i){
            if(i){
                fprintf(out, " ");
            }
            fprintf(out, "%d", cnt+i);
        }
        fprintf(out, "\n");
        cnt+=m;
    }

    fclose(out);
    return true;
}

bool Utils::writeMesh(const std::string &filePath, const Eigen::MatrixXd &V,
                      const Eigen::MatrixXi &C){
    FILE* out = fopen(filePath.c_str(), "w");
    if(!out){
        std::cout << "error on writing mesh" << std::endl;
        return false;
    }

    int nV = V.rows();
    int nEle = C.rows();
    int m = C.cols();
    // write number of C and F
    fprintf(out, "%d\n%d %d\n", nV, nEle, m);

    for(int vI=0; vI<nV; ++vI){
        fprintf(out, "%.18le %.18le %.18le\n", V(vI, 0), V(vI, 1), V(vI, 2));
    }
    for(int eleI=0; eleI<nEle; ++eleI){
        for(int j=0; j < m; ++j){
            if(j){
                fprintf(out, " ");
            }
            fprintf(out, "%d", C(eleI, j));
        }
        fprintf(out, "\n");
    }

    fclose(out);
    return true;
}

} // namespace SIM