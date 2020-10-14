#include "mymath.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <structopt/app.hpp>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
using namespace GiNaC;

using MatrixDn = Matrix<MaxAlgebra, Dynamic, Dynamic>;

MaxAlgebra
spectral_radius(MatrixDn& m)
{
  auto n(m.rows());
  MaxAlgebra radius(m.trace());
  auto m_i(m);

  for (long i = 2; i <= n; ++i) {
    m_i *= m;
    radius += pow(m_i.trace().value, inverse(i));
  }

  return radius;
}

MatrixDn
cleany(MatrixDn& m)
{
  auto n = m.rows();
  auto lambda = spectral_radius(m);
  auto inv_lambda = MaxAlgebra(1) / lambda;

  MatrixDn result = MatrixDn::Identity(n, n);

  MatrixDn A_lambda = inv_lambda * m;
  auto A_i(A_lambda);

  for (long i = 1; i < n; ++i) {
    result += A_i;
    A_i *= A_lambda;
  }

  return result;
}

int
main(int argc, char* argv[])
{

  // Пример решений (стр. 54)

  MatrixDn A(4, 4);
  A << 1, 3, 4, 2,                                     //
    numeric(1) / 3, 1, numeric(1) / 2, numeric(1) / 3, //
    numeric(1) / 4, 2, 1, 4,                           //
    numeric(1) / 2, 3, numeric(1) / 4, 1;              //

  std::cout << A << std::endl;
  MatrixDn cleany_A = cleany(A);
  std::cout << cleany_A << std::endl;

  std::cout << cleany_A.householderQr().matrixQR() << std::endl;

  return 0;
}
