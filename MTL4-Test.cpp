#include <iostream>
#include <boost/numeric/mtl/mtl.hpp>

int main()
{
    using namespace mtl;

    const unsigned n= 10;
    compressed2D<double>                         A(n, n);
    dense2D<int, matrix::parameters<col_major> > B(n, n);
    morton_dense<double, 0x555555f0>             C(n, n), D(n, n);

    matrix::laplacian_setup(A, 2, 5);
    matrix::hessian_setup(B, 1); matrix::hessian_setup(C, 2.0); matrix::hessian_setup(D, 3.0);

    D+= A - 2 * B + C;

    std::cout << "The matrices are: A=\n" << A << "B=\n" << B << "C=\n" << C << "D=\n" << D;

    return 0;
}
