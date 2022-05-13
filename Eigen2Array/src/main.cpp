#include <iostream>
#include <Eigen/Eigen>

void eigenVec2array(){
    Eigen::VectorXd e(3);
    e << 1, 2, 3;
    int len = e.size();
    double *a = e.data(); // e and a point to the same memory
    a[0] = -1;

    std::cout << "e\n" << e << std::endl << std::endl;
    std::cout << "a\n";
    for(int i=0; i<len; ++i){
        std::cout << a[i] << std::endl;
    }

    double *a2 = new double[len];
    Eigen::Map<Eigen::VectorXd>(a2, len) = e; // e and a2 point to different memory
    a2[0] = -2;

    std::cout << "\ne\n" << e << std::endl << std::endl;
    std::cout << "a2\n";
    for(int i=0; i<len; ++i){
        std::cout << a2[i] << std::endl;
    }
}

void eigenMat2array(){
    Eigen::MatrixXd e(3, 3);
    e << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    int row = e.rows(), col = e.cols();
    double *a = e.data(); // e and a point to the same memory
    a[0] = -1;

    std::cout << "e\n" << e << std::endl << std::endl;
    std::cout << "a\n";
    for(int i=0; i<row*col; ++i){
        std::cout << a[i] << std::endl;
    }

    double *a2 = new double[row*col];
    Eigen::Map<Eigen::MatrixXd>(a2, row, col) = e; // e and a2 point to different memory
    a2[0] = -2;

    std::cout << "e\n" << e << std::endl << std::endl;
    std::cout << "a2\n";
    for(int i=0; i<row*col; ++i){
        std::cout << a2[i] << std::endl;
    }
}

void array2eigenVec(){
    int len = 3;
    double a[len];
    for(int i=0; i<len; ++i){
        a[i] = i + 1;
    }
    Eigen::Map<Eigen::VectorXd> e(a, len); // e and a point to the same memory
    e(0) = -1;

    std::cout << "a\n";
    for(int i=0; i<len; ++i){
        std::cout << a[i] << std::endl;
    }
    std::cout << "e\n" << e << std::endl;

    Eigen::VectorXd e2 = Eigen::Map<Eigen::VectorXd>(a, len); // e2 and a point to different memory
    e2(0) = -2;

    std::cout << "\na\n";
    for(int i=0; i<len; ++i){
        std::cout << a[i] << std::endl;
    }
    std::cout << "e2\n" << e2 << std::endl;
}

void array2eigenMat(){
    int row = 3, col = 3;
    double a[row*col];
    a[0] = 1, a[3] = 2, a[6] = 3;
    a[1] = 4, a[4] = 5, a[7] = 6;
    a[2] = 7, a[5] = 8, a[8] = 9;
    Eigen::Map<Eigen::MatrixXd> e(a, row, col); // e and a point to the same memory
    e(0, 0) = -1;

    std::cout << "e\n" << e << std::endl << std::endl;
    std::cout << "a\n";
    for(int i=0; i<row*col; ++i){
        std::cout << a[i] << std::endl;
    }

    // e2 and a point to different memory
    Eigen::MatrixXd e2 = Eigen::Map<Eigen::MatrixXd>(a, row, col);
    e(0, 0) = -2;

    std::cout << "\ne2\n" << e2 << std::endl << std::endl;
    std::cout << "a\n";
    for(int i=0; i<row*col; ++i){
        std::cout << a[i] << std::endl;
    }
}

int main(int argc, char** argv){
    // eigenVec2array();
    // eigenMat2array();
    // array2eigenVec();
    array2eigenMat();

    return 0;
}