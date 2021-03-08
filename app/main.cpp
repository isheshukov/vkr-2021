#include <Eigen/Dense>
#include <iostream>
#include <structopt/app.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include "mymath.hpp"

using MatrixDn = Eigen::Matrix<MaxAlgebra, Eigen::Dynamic, Eigen::Dynamic>;
using ma = MaxAlgebra;

static Eigen::IOFormat StdoutFmt(4, 0, ", ", "\n", "\t", "", "\n", "");
static Eigen::IOFormat LaTeXFmt(Eigen::FullPrecision, Eigen::DontAlignCols,
                                "& ", "\\\\", "", "", "\\begin{pmatrix}",
                                "\\end{pmatrix}");
static Eigen::IOFormat MatlabFmt(Eigen::FullPrecision, Eigen::DontAlignCols,
                                 ", ", "\n ", "[", "]", "", "");
static Eigen::IOFormat MathematicaFmt(Eigen::FullPrecision,
                                      Eigen::DontAlignCols, ", ", ",\n ", "{",
                                      "}", "{", "}");

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
    // std::cout << "\nm^" << i << "= " << m_i.format(LaTeXFmt) << "\n";
  }

  return radius;
}
MatrixDn cleany(MatrixDn& m) {
  // std::cout << "\n\tCLEANY MATRIX";
  // std::cout << "\n\t---------------";

  auto n = m.rows();

  MatrixDn result = MatrixDn::Identity(n, n);
  // std::cout << "\n\tI=\n" << result.format(MathematicaFmt) << "\n";
  // std::cout << "\n\tI= " << result.format(LaTeXFmt) << "\n";

  auto A_i(m);

  for (long i = 1; i < n; ++i) {
    result += A_i;
    // std::cout << "\n\tm^" << i << "=\n" << A_i.format(LaTeXFmt) << "\n";
    A_i *= m;
  }

  return result;
}
std::pair<MatrixDn::Index, MatrixDn::Index> best_diff_vector_coefficients(
    MatrixDn& P) {
  std::vector<MaxAlgebra> k_candidates;
  auto rows = P.rows();
  auto cols = P.cols();
  auto ones_vec = MatrixDn::Ones(rows, 1);

  for (size_t i = 0; i < cols; ++i) {
    auto p = P.col(i);
    auto pminus = p.transpose().cwiseInverse();
    auto pre_result = ones_vec.transpose() * p * pminus * ones_vec;
    std::cout << pre_result << std::endl;
    k_candidates.push_back(pre_result.eval()(0, 0));
  }

  auto k =
      std::distance(k_candidates.begin(),
                    std::max_element(k_candidates.begin(), k_candidates.end()));

  auto p_k_inv = P.col(k).cwiseInverse().eval();
  MatrixDn::Index l_index;
  p_k_inv.maxCoeff(&l_index);

  return {k, l_index};
}

int main(int argc, char* argv[]) {
  // Пример решений (стр. 54)
  MatrixDn C(6, 6);
  C << 1, MaxAlgebra(1, 4), MaxAlgebra(1, 5), MaxAlgebra(1, 4), 6,
      MaxAlgebra(1, 6), 4, 1, MaxAlgebra(1, 3), 3, 6, MaxAlgebra(1, 2), 5, 3, 1,
      4, 7, 3, 4, MaxAlgebra(1, 3), MaxAlgebra(1, 4), 1, 5, MaxAlgebra(1, 5),
      MaxAlgebra(1, 6), MaxAlgebra(1, 6), MaxAlgebra(1, 7), MaxAlgebra(1, 5), 1,
      MaxAlgebra(1, 7), 6, 2, MaxAlgebra(1, 3), 5, 7, 1;

  MatrixDn A1(3, 3);
  MatrixDn A2(3, 3);
  MatrixDn A3(3, 3);
  MatrixDn A4(3, 3);
  MatrixDn A5(3, 3);
  MatrixDn A6(3, 3);

  auto lambda = spectral_radius(C);

  MatrixDn Calmoststar = GiNaC::ex(1) / GiNaC::ex(lambda) * C;
  auto w = cleany(Calmoststar);
  std::cout << "w = " << w.format(MathematicaFmt) << std::endl;

  auto delta = MatrixDn::Ones(1, 6) * w * MatrixDn::Ones(6, 1);
  std::cout << "delta = " << delta << std::endl;
  std::cout << "delta^-1 = " << 1 / delta.value() << std::endl;

  MatrixDn W1m = (1 / delta.value() * MatrixDn::Ones(6, 6)) + Calmoststar;
  auto w1 = cleany(W1m);

  MatrixDn P(6, 2);
  P << w.col(0), w.col(1);
  std::cout << "P = " << P.format(MathematicaFmt) << std::endl;

  auto [k, l] = best_diff_vector_coefficients(P);

  MatrixDn Plk_inverse = MatrixDn::Zero(P.rows(), P.cols());
  Plk_inverse(l, k) = 1 / P(l, k);

  auto w2 = P * (MatrixDn::Identity(P.cols(), P.cols()) +
                 Plk_inverse.transpose() * P);

  std::cout << "w2 = " << w2.format(MathematicaFmt) << std::endl;
  return 0;
}
