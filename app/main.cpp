#include <Eigen/Dense>
#include <iostream>
#include <structopt/app.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include "mymath.hpp"

using MatrixDn = Eigen::Matrix<MaxAlgebra, Eigen::Dynamic, Eigen::Dynamic>;

MaxAlgebra spectral_radius(MatrixDn& m) {
  auto n(m.rows());
  MaxAlgebra radius(m.trace());
  auto m_i(m);

  //  std::cout << "\nSPECTRAL RADIUS";
  //  std::cout << "\n---------------";

  //  std::cout << "\nm^" << 1 << "= " << m_i << "\n";

  for (long i = 2; i <= n; ++i) {
    m_i *= m;
    radius += pow(m_i.trace().value, GiNaC::inverse(i));
    //    std::cout << "\nm^" << i << "= " << m_i << "\n";
  }

  return radius;
}

MatrixDn cleany(MatrixDn& m) {
  // std::cout << "\nCLEANY MATRIX";
  // std::cout << "\n---------------";

  auto n = m.rows();
  auto lambda = spectral_radius(m);
  auto inv_lambda = MaxAlgebra(1) / lambda;

  MatrixDn result = MatrixDn::Identity(n, n);
  // std::cout << "\nm^" << 0 << "=\n" << result << "\n";

  MatrixDn A_lambda = inv_lambda * m;
  auto A_i(A_lambda);

  for (long i = 1; i < n; ++i) {
    result += A_i;
    A_i *= A_lambda;
    // std::cout << "\nm^" << i << "=\n" << result << "\n";
  }

  return result;
}

int main(int argc, char* argv[]) {
  // Пример решений (стр. 54)
  Eigen::IOFormat StdoutFmt(4, 0, ", ", "\n", "\t", "", "\n", "");
  Eigen::IOFormat LaTeXFmt(Eigen::FullPrecision, Eigen::DontAlignCols, "& ",
                           "\\\\", "", "", "\\begin{pmatrix}",
                           "\\end{pmatrix}");

  MatrixDn C(2, 2);
  MatrixDn A1(3, 3);
  MatrixDn A2(3, 3);
  C << 1, MaxAlgebra(1, 3), 3, 1;
  A1 << 1, 3, GiNaC::ex(1) / 3, GiNaC::ex(1) / 3, 1, 1, 3, 1, 1;
  A2 << 1, GiNaC::ex(1) / 3, 5, 3, 1, 7, GiNaC::ex(1) / 5, GiNaC::ex(1) / 7, 1;

  auto lambda = spectral_radius(C);
  std::cout << "Spectral radius of C (lambda) = " << lambda << std::endl;

  auto B = (A1 + 3 * A2).eval();
  std::cout << "Matrix B is " << B.format(StdoutFmt) << std::endl;

  auto mu = spectral_radius(B);
  std::cout << "Spectral radius of B (mu) = " << mu << std::endl;

  auto delta = GiNaC::ex(7) / GiNaC::sqrt(GiNaC::ex(5));

  auto D = MatrixDn::Ones(3, 3) / GiNaC::ex(delta);
  auto BB = (1 / mu * B).eval();
  auto BBstar = cleany(BB);
  std::cout << "Cleany matrix of 1 / mu * B\n"
            << BBstar.format(StdoutFmt) << std::endl;

  auto M = (D + BB).eval();
  std::cout << "Matrix (D + 1/mu * B) is " << M.format(StdoutFmt) << std::endl;
  std::cout << "Cleany matrix of D + 1 / mu * B" << cleany(M).format(StdoutFmt)
            << std::endl;

  return 0;
}
