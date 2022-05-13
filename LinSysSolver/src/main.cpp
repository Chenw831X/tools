#include <iostream>
#include <chrono>
#include <Eigen/Eigen>

#ifdef LINSYSSOLVER_USE_CHOLMOD
#include "CHOLMODSolver.hpp"
#else
#include "EigenLibSolver.hpp"
#endif

int N;
int nnz;

void readF(const std::string &filePath, Eigen::VectorXd &rhs){
    FILE* in = fopen(filePath.c_str(), "r");
    assert(in);

    int numRows = N;
    // std::cout << "number of numRows: " << numRows << std::endl;
    rhs.resize(numRows);
    for(int i=0; i<numRows; ++i){
        fscanf(in, "%le", &rhs[i]);
    }

    fclose(in);
    std::cout << "read F done" << std::endl;
}

void readMatrix(const std::string &filePath, Eigen::SparseMatrix<double> &Sp){
    FILE* in = fopen(filePath.c_str(), "r");
    assert(in);

    Eigen::VectorXi I, J;
    Eigen::VectorXd A;

    I.resize(nnz);
    for(int i=0; i<nnz; ++i){
        int tmp;
        fscanf(in, "%d", &tmp);
        I(i) = tmp - 1;
    }

    J.resize(nnz);
    for(int i=0; i<nnz; ++i){
        int tmp;
        fscanf(in, "%d", &tmp);
        J(i) = tmp - 1;
    }

    A.resize(nnz);
    for(int i=0; i<nnz; ++i){
        fscanf(in, "%le", &A[i]);
    }

    std::vector<Eigen::Triplet<double>> triLists;
    triLists.reserve(nnz);
    for(int i=0; i<nnz; ++i){
        triLists.emplace_back(Eigen::Triplet<double>(I(i), J(i), A(i)));
    }
    Sp.setFromTriplets(triLists.begin(), triLists.end());
    // std::cout << Sp << std::endl;

    fclose(in);
    std::cout << "read K done" << std::endl;
}

void loadF(const char* filePath, Eigen::VectorXd &rhs){
    FILE* in = fopen(filePath, "rb");
    assert(in);

    fread(&N, sizeof(N), 1, in);
    rhs.resize(N);
    fread(rhs.data(), sizeof(rhs[0]), N, in);

    fclose(in);
    std::cout << "load F done" << std::endl;
}

void loadK(const char* filePath, Eigen::SparseMatrix<double> &Sp){
    FILE* in = fopen(filePath, "rb");
    assert(in);

    fread(&nnz, sizeof(nnz), 1, in);
    nnz /= 3;
    Eigen::VectorXi ia, ja;
    Eigen::VectorXd a;

    ia.resize(nnz);
    fread(ia.data(), sizeof(ia[0]), nnz, in);

    ja.resize(nnz);
    fread(ja.data(), sizeof(ja[0]), nnz, in);

    ia.array() -= 1;
    ja.array() -= 1;

    a.resize(nnz);
    fread(a.data(), sizeof(a[0]), nnz, in);

    std::vector<Eigen::Triplet<double>> lists;
    lists.reserve(nnz);
    for(int i=0; i < nnz; ++i){
        lists.emplace_back(ia(i), ja(i), a(i));
    }
    Sp.setFromTriplets(lists.begin(), lists.end());

    fclose(in);
    std::cout << "load K done" << std::endl;
}

bool writeMatrix(const std::string &filePath, const Eigen::MatrixXd &mat){
    FILE* out = fopen(filePath.c_str(), "w");
    if(!out){
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

void writeK(Eigen::SparseMatrix<double> &Sp)
{
    FILE* outI = fopen("../output/I.txt", "w");
    assert(outI);
    FILE* outJ = fopen("../output/J.txt", "w");
    assert(outJ);
    FILE* outA = fopen("../output/A.txt", "w");
    assert(outA);

    for (int k=0; k<Sp.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Sp, k); it; ++it) {
            fprintf(outI, "%d\n", (int)(it.row()) + 1);
            fprintf(outJ, "%d\n", (int)(it.col()) + 1);
            fprintf(outA, "%.18le\n", it.value());
            // it.value();
            // it.row();   // row index
            // it.col();   // col index (here it is equal to k)
            // it.index(); // inner index, here it is equal to it.row()
        }
    }

    fclose(outI);
    fclose(outJ);
    fclose(outA);
}

int main(){
    Eigen::VectorXd rhs;
    loadF("../input/F_.bin", rhs);
    std::cout << "N: " << N << std::endl;
    Eigen::SparseMatrix<double> Sp(N, N);
    loadK("../input/K_.bin", Sp);

    std::cout << "nnz: " << Sp.nonZeros() << std::endl;

    SIM::LinSysSolver<Eigen::VectorXi, Eigen::VectorXd> *linSysSolver;

#ifdef LINSYSSOLVER_USE_CHOLMOD
    linSysSolver = new SIM::CHOLMODSolver<Eigen::VectorXi, Eigen::VectorXd>();
#else
    linSysSolver = new SIM::EigenLibSolver<Eigen::VectorXi, Eigen::VectorXd>();
#endif

    linSysSolver->set_pattern(Sp);
    linSysSolver->analyze_pattern();

    auto start = std::chrono::steady_clock::now();

    linSysSolver->factorize();
    Eigen::VectorXd solution;
    linSysSolver->solve(rhs, solution);

    std::chrono::duration<double> time_(std::chrono::steady_clock::now() - start);
    std::cout << "time: " << time_.count() << std::endl;

    Eigen::VectorXd KU;
    linSysSolver->multiply(solution, KU);
    std::cout << (KU - rhs).cwiseAbs().mean() << std::endl;
    std::cout << (KU - rhs).cwiseAbs().maxCoeff() << std::endl;
    // std::cout << "U.mean: " << solution.mean() << std::endl;
    std::cout << "U * F: " << solution.transpose() * rhs << std::endl;
    // Eigen::VectorXd U_m;
    // readF("../output/U_SIZE_179046x1.txt", U_m);
    // std::cout << "U_M * F:" << U_m.transpose() * rhs << std::endl;
    // std::cout << "U_m.mean: " << U_m.mean() << std::endl;

    writeMatrix("../output/U_cho.txt", solution);


//     SIM::LinSysSolver<Eigen::VectorXi, Eigen::VectorXd> *linSysSolver;
//
// #ifdef LINSYSSOLVER_USE_CHOLMOD
//     linSysSolver = new SIM::CHOLMODSolver<Eigen::VectorXi, Eigen::VectorXd>();
// #else
//     linSysSolver = new SIM::EigenLibSolver<Eigen::VectorXi, Eigen::VectorXd>();
// #endif
//
//     Eigen::VectorXd rhs, solution;
//     linSysSolver->load("../input/K.txt", rhs);
//     linSysSolver->analyze_pattern();
//
//     auto start = std::chrono::steady_clock::now();
//
//     linSysSolver->factorize();
//     linSysSolver->solve(rhs, solution);
//
//     std::chrono::duration<double> time_(std::chrono::steady_clock::now() - start);
//     std::cout << "time: " << time_.count() << std::endl;
//
//     writeMatrix("../output/U_cho.txt", solution);

}