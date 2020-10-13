#include "mymath.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <structopt/app.hpp>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
using namespace GiNaC;

using Matrix44 = Matrix<MaxAlgebra, 4, 4>;
using MatrixDn = Matrix<MaxAlgebra, Dynamic, Dynamic>;
using M = MaxAlgebra;

numeric
spectral_radius(MatrixDn& m)
{
  auto n = m.rows();
  auto radius(m.trace().value);
  auto m_i = m;
  for (long i = 2; i <= n; ++i) {
    m_i *= m;
    auto mitr = m_i.trace().value;
    auto tr = ex_to<numeric>(pow(mitr, numeric(1, i)).eval());
    radius = std::max(radius, tr);
  }

  std::cout << "radius " << radius << std::endl;
  return radius;
}

MatrixDn
cleany(MatrixDn& m)
{
  auto n = m.rows();
  auto lambda = ex_to<numeric>(spectral_radius(m).eval());
  MaxAlgebra inv_lambda(ex_to<numeric>(1 / lambda));

  MatrixDn result(n, n);
  for (long i = 0; i < n; ++i) {
    for (long j = 0; j < n; ++j) {
      result(i, j) = MaxAlgebra(0);
    }
  }

  for (long i = 0; i < n; ++i) {
    result(i, i) = MaxAlgebra(1);
  }

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
  A << M(1), M(3), M(4), M(2), M(1, 3), M(1), M(1, 2), M(1, 3), M(1, 4), M(2),
    M(1), M(4), M(1, 2), M(3), M(1, 4), M(1);

  std::cout << cleany(A) << std::endl;

  return 0;
}
