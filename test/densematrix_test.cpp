#include "../src/types/include/densematrix.hpp"

#include <iostream>

int main()
{
    std::cout << "FMa" << std::endl;
    hptypes::DenseMatrix FMa(4,4);
    FMa.set_local(0, 0, 2);
    FMa.set_local(0, 1,-3);
    FMa.set_local(0, 2, 0);
    FMa.set_local(0, 3, 3);
    FMa.set_local(1, 0, 2);
    FMa.set_local(1, 1, 4);
    FMa.set_local(1, 2, 3);
    FMa.set_local(1, 3, 1);
    FMa.set_local(2, 0,-3);
    FMa.set_local(2, 1, 6);
    FMa.set_local(2, 2,-3);
    FMa.set_local(2, 3, 0);
    FMa.set_local(3, 0, 2);
    FMa.set_local(3, 1, 2);
    FMa.set_local(3, 2, 3);
    FMa.set_local(3, 3, 4);
    FMa.print_local();
    std::cout << std::endl;
    std::cout << "FMainv" << std::endl;
    hptypes::DenseMatrix FMainv(FMa.get_inverse());
    FMainv.print_local();
    std::cout << std::endl;
    std::cout << "FMa * FMainv" << std::endl;
    FMa.get_mat_mul(FMainv).print_local();

//    FMa.print_local();
//    std::cout << std::endl;
//    std::cout << "FMb" << std::endl;
//    hptypes::DenseMatrix<double> FMb(FMa.get_transpose());
//    FMb.print_local();
//    std::cout << std::endl;
//    FMb = FMb.get_transpose();
//    FMb.print_local();
//    hptypes::DenseMatrix<double> FMb(4,4);
//    FMb.set_local(0, 0, 1);
//    FMb.set_local(2, 0, 2);
//    FMb.set_local(0, 2, 2);
//    FMb.set_local(0, 3, 3);
//    std::cout << FMb.issymmetric() << std::endl;
    return 0;
}
