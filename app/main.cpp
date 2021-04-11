#include <Eigen/Dense>
#include <iostream>
#include <structopt/app.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include "mymath.hpp"

/*
 * MaxAlgebra -- класс, реализующий арифметику в max-алгебре.
 * MatrixDn -- матрица, скаляры которой находятся в max-алгебре.
 * ma -- просто сокращение для удобства написания.
 */
using MatrixDn = Eigen::Matrix<MaxAlgebra, Eigen::Dynamic, Eigen::Dynamic>;
using ma = MaxAlgebra;

/*
 * Данные строки отвечают за различные способы вывода матриц и векторов.
 *      StdoutFmt -- обычное отображение элементов через пробел
 *      LaTeXFmt -- формат, пригодный для вставки в tex-файл.
 *      MatlabFmt -- формат, пригодный для вставки в matlab-скрипт.
 *      MathematicaFmt -- формат, пригодный для вставки в mathematica-скрипт.
 *
 * Использование:
 *
 *      std::cout << w.format(format) << std::endl;
 *
 * где w -- матрица.
 */
static Eigen::IOFormat StdoutFmt(4, 0, ", ", "\n", "\t", "", "\n", "");
static Eigen::IOFormat LaTeXFmt(Eigen::FullPrecision, Eigen::DontAlignCols,
                                "& ", "\\\\", "", "", "\\begin{pmatrix}",
                                "\\end{pmatrix}");
static Eigen::IOFormat MatlabFmt(Eigen::FullPrecision, Eigen::DontAlignCols,
                                 ", ", "\n ", "[", "]", "", "");
static Eigen::IOFormat MathematicaFmt(Eigen::FullPrecision,
                                      Eigen::DontAlignCols, ", ", ",\n ", "{",
                                      "}", "{", "}");

MatrixDn matrixWithoutColumn(MatrixDn& matrix, unsigned int colToRemove)
{
    auto result = matrix;
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        result.block(0,colToRemove,numRows,numCols-colToRemove) = result.block(0,colToRemove+1,numRows,numCols-colToRemove);
    result.conservativeResize(numRows,numCols);

    return result.eval();
}

/*
 * Функция вычисления спектрального радиуса матрицы.
 */
MaxAlgebra spectral_radius(MatrixDn& m) {
  auto n(m.rows());
  MaxAlgebra radius(m.trace());
  auto m_i(m);

  for (long i = 2; i <= n; ++i) {
    m_i *= m;
    radius += pow(m_i.trace().value, GiNaC::inverse(i));
  }

  return radius;
}

MatrixDn pseudo_inverse(const MatrixDn& A)
{
    return A.transpose().cwiseInverse().eval();
}

/*
 * Функция невязки.
 * Delta(A, b) = ...
 */
MaxAlgebra residual_function(MatrixDn& A, Eigen::Matrix<MaxAlgebra, -1, 1, 0, -1, 1>& b) {
    auto Api = pseudo_inverse(A);
    auto bpi = pseudo_inverse(b);

    return (pseudo_inverse(A*pseudo_inverse((bpi*A).eval()))*b).value();
}

MatrixDn linearly_independent_system(MatrixDn& A){
    auto result = A;

    for (int i = 0; i < A.cols(); ++i){
        auto matrix_without_i_col = matrixWithoutColumn(A, i);
        auto col = A.col(i).eval();
        if (residual_function(matrix_without_i_col, col) == 1) {
            // delete vector
            result = matrix_without_i_col;
        }
    }

    return result;
}

/*
 * Функция, вычисления матрицы Клини.
 */
MatrixDn cleany(MatrixDn& m) {
  auto n = m.rows();

  MatrixDn result = MatrixDn::Identity(n, n);
  auto A_i(m);

  for (long i = 1; i < n; ++i) {
    result += A_i;
    A_i *= m;
  }

  return result;
}

/*
 * Функция для получения коэфициентов k и l в задаче нахождения наилучшего
 * дифференцирующего вектора.
 */
std::pair<MatrixDn::Index, MatrixDn::Index> best_diff_vector_coefficients(
    MatrixDn& P) {
  std::vector<MaxAlgebra> k_candidates;
  auto rows = P.rows();
  auto cols = P.cols();
  auto ones_vec = MatrixDn::Ones(rows, 1);

  for (long i = 0; i < cols; ++i) {
    auto p = P.col(i);
    auto pminus = p.transpose().cwiseInverse();
    auto pre_result = ones_vec.transpose() * p * pminus * ones_vec;
    k_candidates.push_back(pre_result.value());
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
  /*
   * Задание исходных данных.
   *
   * Примеры:
   *    Создание матрицы 2x2 с целочисленными элементами.
   *        MatrixDn A(2,2);
   *        A << 1, 2, 3, 4;
   *
   *    Создание матрицы 2x2 с дробными элементами элементами.
   *        MatrixDn A(2,2);
   *        A << MaxAlgebra(1, 2), 2, MaxAlgebra(1, 3), 4;
   *
   *    Запись 3^2
   *        pow(MaxAlgebra(3), 2);
   *
   *    1/(3^2)
   *        1 / pow(MaxAlgebra(3), 2))
   *
   *    Матрица из единиц, нулей размерности N
   *        MatrixDn::Ones(N, N);
   *        MatrixDn::Zero(N, N);
   *
   *    Единичная матрица размерности N
   *        MatrixDn::Identity(N, N);
   */
  MatrixDn C(6, 6);
  C << 1, MaxAlgebra(1, 4), MaxAlgebra(1, 5), MaxAlgebra(1, 4), 6,
      MaxAlgebra(1, 6), 4, 1, MaxAlgebra(1, 3), 3, 6, MaxAlgebra(1, 2), 5, 3, 1,
      4, 7, 3, 4, MaxAlgebra(1, 3), MaxAlgebra(1, 4), 1, 5, MaxAlgebra(1, 5),
      MaxAlgebra(1, 6), MaxAlgebra(1, 6), MaxAlgebra(1, 7), MaxAlgebra(1, 5), 1,
      MaxAlgebra(1, 7), 6, 2, MaxAlgebra(1, 3), 5, 7, 1;

  MatrixDn A1(3, 3);
  A1 << 1, 5, 8, ma(1, 5), 1, 5, ma(1, 8), ma(1, 5), 1;
  MatrixDn A2(3, 3);
  A2 << 1, 7, 9, ma(1, 7), 1, 7, ma(1, 9), ma(1, 7), 1;
  MatrixDn A3(3, 3);
  A3 << 1, ma(1, 7), ma(1, 9), 7, 1, ma(1, 7), 9, 7, 1;
  MatrixDn A4(3, 3);
  A4 << 1, 3, 5, ma(1, 3), 1, 4, ma(1, 5), ma(1, 4), 1;
  MatrixDn A5(3, 3);
  A5 << 1, 3, 5, ma(1, 3), 1, 4, ma(1, 5), ma(1, 4), 1;
  MatrixDn A6(3, 3);
  A6 << 1, 7, 9, ma(1, 7), 1, 7, ma(1, 9), ma(1, 7), 1;

  /*
   * Вычисление спектрального радиуса
   */
  std::cout << "Вычисление спектрального радиуса" << std::endl;
  auto lambda = spectral_radius(C);
  std::cout << "lambda = " << lambda << std::endl << std::endl;

  /*
   * Вычисление (\lambda^{-1}C)^*
   */
  std::cout << "Вычисление (\\lambda^{-1}C)^*" << std::endl;
  MatrixDn Calmoststar = 1 / lambda * C;
  auto w = cleany(Calmoststar);
  std::cout << "w = " << w.format(MatlabFmt) << std::endl << std::endl;

  /*
   * Вычисление \delta
   */
  std::cout << "Вычисление \\delta" << std::endl;
  auto delta = (MatrixDn::Ones(1, 6) * w * MatrixDn::Ones(6, 1)).value();
  // .value() нужно, чтобы достать скаляр из матрицы 1x1.
  std::cout << "delta = " << delta << std::endl;
  std::cout << "delta^-1 = " << 1 / delta << std::endl << std::endl;

  /*
   * Вычисление наихудшего дифференцирующего вектора
   */
  std::cout << "Вычисление наихудшего дифференцирующего вектора" << std::endl;
  MatrixDn W1m = (1 / delta * MatrixDn::Ones(6, 6)) + Calmoststar;
  auto w1 = cleany(W1m);
  std::cout << "w1 = " << w1.format(MatlabFmt) << std::endl << std::endl;

  /* Получилась матрица из двух ЛН векторов. Попробуем получить
   * "самый наихудший" вектор.
   * Возьмём первые два столбца (я знаю, что они ЛН), отнормируем по максимуму
   * и сложим.
   */

  auto worst_vector =
      (w.col(0) / w.col(0).maxCoeff()) + (w.col(1) / w.col(1).maxCoeff());
  std::cout << "worst_vector = " << worst_vector.format(MatlabFmt) << std::endl
            << std::endl;

  /*
   * Получим взвешенную сумму.
   */

  auto B1 =
      (worst_vector(0) * A1 + worst_vector(1) * A2 + worst_vector(2) * A3 +
       worst_vector(3) * A4 + worst_vector(4) * A5 + worst_vector(5) * A6)
          .eval();

  auto nu = spectral_radius(B1);
  std::cout << "B1 = " << B1.format(MatlabFmt) << std::endl;
  std::cout << "nu = " << nu << std::endl;

  MatrixDn almost_x1 = 1 / nu * B1;
  auto x1 = cleany(almost_x1);

  std::cout << "x1 = " << x1.format(MatlabFmt) << std::endl;

  /*
   * Случай неединственности решения.
   */

  auto delta1 = (MatrixDn::Ones(1, 3) * x1 * MatrixDn::Ones(3, 1)).value();
  MatrixDn almost_x1new = (delta1 * MatrixDn::Ones(3, 3)) + almost_x1;
  auto x1new = cleany(almost_x1new);
  std::cout << "x1new = " << x1new.format(MatlabFmt) << std::endl;

  /*
   * Вычисление наилучшего дифференцирующего вектора.
   * Матрица P --- линейно-независимые столбцы матрицы w.
   * Временно, получение линейно-независимых столбцов происходит
   * при помощи отдельного скрипта Wolfram Mathematica.
   */
  std::cout << "Вычисление наилучшего дифференцирующего вектора" << std::endl;
  MatrixDn P(6, 2);
  P << w.col(0), w.col(1);
  std::cout << "P = " << P.format(MatlabFmt) << std::endl << std::endl;

  /*
   * Получение коэфициентов k и l.
   */
  std::cout << "Получение коэфициентов k и l" << std::endl;
  auto [k, l] = best_diff_vector_coefficients(P);
  std::cout << "k = " << k << ", l = " << l << std::endl << std::endl;

  /*
   * Вычисление w_2.
   */
  std::cout << "Вычисление w_2" << std::endl;
  MatrixDn Plk_inverse = MatrixDn::Zero(P.rows(), P.cols());
  Plk_inverse(l, k) = 1 / P(l, k);

  auto w2 = P * (MatrixDn::Identity(P.cols(), P.cols()) +
                 Plk_inverse.transpose() * P);

  std::cout << "w2 = " << w2.format(MatlabFmt) << std::endl;
  return 0;
}
