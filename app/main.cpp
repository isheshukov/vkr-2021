#include <Eigen/Dense>
#include <iostream>
#include <structopt/app.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include "mymath.hpp"

/*
 * MaxTimes -- класс, реализующий арифметику в max-алгебре.
 * MatrixDn -- матрица, скаляры которой находятся в max-алгебре.
 * ma -- просто сокращение для удобства написания.
 */
using MatrixDn = Eigen::Matrix<MaxTimes, Eigen::Dynamic, Eigen::Dynamic>;
using ma = MaxTimes;

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
MaxTimes spectral_radius(MatrixDn& m) {
  auto n(m.rows());
  MaxTimes radius(m.trace());
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
MaxTimes residual_function(MatrixDn& A, Eigen::Matrix<MaxTimes, -1, 1, 0, -1, 1>& b) {
    auto Api = pseudo_inverse(A);
    auto bpi = pseudo_inverse(b);

    return (pseudo_inverse(A*pseudo_inverse((bpi*A).eval()))*b).value();
}

MatrixDn linearly_independent_system(MatrixDn& A){
    auto result = A;
    auto num_cols = A.cols();
    size_t cur_col = 0;

    for (int i = 0; i < num_cols; ++i){
        auto matrix_without_i_col = matrixWithoutColumn(result, cur_col);
        auto col = result.col(cur_col).eval();
        if (residual_function(matrix_without_i_col, col) == 1) {
            result = matrix_without_i_col;
        } else {
            cur_col++;
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
std::vector<std::pair<MatrixDn::Index, MatrixDn::Index>> best_diff_vector_coefficients(
    MatrixDn& P) {
  std::vector<std::pair<MatrixDn::Index, MatrixDn::Index>> result;
  std::vector<MaxTimes> k_candidates;
  auto rows = P.rows();
  auto cols = P.cols();
  auto ones_vec = MatrixDn::Ones(rows, 1);

  for (long i = 0; i < cols; ++i) {
    auto p = P.col(i);
    auto pminus = p.transpose().cwiseInverse();
    auto pre_result = ones_vec.transpose() * p * pminus * ones_vec;
    k_candidates.push_back(pre_result.value());
  }

  auto max_k_value = *std::max_element(k_candidates.begin(), k_candidates.end());
  for (MatrixDn::Index k = 0; k < k_candidates.size(); ++k){
      if (max_k_value == k_candidates[k]) {
          auto p_k_inv = P.col(k).cwiseInverse().eval();
          auto max_l_value = p_k_inv.maxCoeff();
          for (MatrixDn::Index l = 0; l < p_k_inv.size(); ++l){
              if (max_l_value == p_k_inv(l)){
                  result.push_back({k, l});
              }
          }
      }
  }


  return result;
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
   *        A << MaxTimes(1, 2), 2, MaxTimes(1, 3), 4;
   *
   *    Запись 3^2
   *        pow(MaxTimes(3), 2);
   *
   *    1/(3^2)
   *        1 / pow(MaxTimes(3), 2))
   *
   *    Матрица из единиц, нулей размерности N
   *        MatrixDn::Ones(N, N);
   *        MatrixDn::Zero(N, N);
   *
   *    Единичная матрица размерности N
   *        MatrixDn::Identity(N, N);
   */

  // Попросить данные
  MatrixDn C(6, 6);
  C << 1, MaxTimes(1, 4), MaxTimes(1, 5), MaxTimes(1, 4), 6,
      MaxTimes(1, 6), 4, 1, MaxTimes(1, 3), 3, 6, MaxTimes(1, 2), 5, 3, 1,
      4, 7, 3, 4, MaxTimes(1, 3), MaxTimes(1, 4), 1, 5, MaxTimes(1, 5),
      MaxTimes(1, 6), MaxTimes(1, 6), MaxTimes(1, 7), MaxTimes(1, 5), 1,
      MaxTimes(1, 7), 6, 2, MaxTimes(1, 3), 5, 7, 1;


  std::vector<MatrixDn>As;
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

  As.push_back(A1);
  As.push_back(A2);
  As.push_back(A3);
  As.push_back(A4);
  As.push_back(A5);
  As.push_back(A6);

//  Маленькая пример из текста диплома
//  MatrixDn C(2,2);
//  C << 1, ma(1,3), 3, 1;

//  MatrixDn A1(3,3);
//  A1 << 1, 3, ma(1, 3),
//        ma(1,3), 1, 1,
//          3, 1, 1;
//  MatrixDn A2(3,3);
//  A2 << 1, 1/3, 5,
//        3, 1, 7,
//        ma(1,5), ma(1,7), 1;

//  As.push_back(A1);
//  As.push_back(A2);

  const auto num_comp_mat = A1.rows();
  const auto num_crit_mat = C.rows();

  /*
   * Вычисление спектрального радиуса
   */
  std::cout << "Computing spectrail radius: ";
  auto lambda = spectral_radius(C);
  std::cout << "lambda = " << lambda << std::endl;

  /*
   * Вычисление (\lambda^{-1}C)^*
   */
  std::cout << "Computing (1/lambda * C)^*" << std::endl;
  MatrixDn Calmoststar = (1 / lambda) * C;
  auto w = cleany(Calmoststar);
  std::cout << "(1/lambda * C)^* = \n" << w.format(MatlabFmt) << std::endl << std::endl;

  /*
   * Вычисление \delta
   */
  std::cout << "Computing delta: " << std::endl;
  auto delta = (MatrixDn::Ones(1, num_crit_mat) * w * MatrixDn::Ones(num_crit_mat, 1)).value();
  // .value() нужно, чтобы достать скаляр из матрицы 1x1.
  std::cout << "delta = " << delta;
  std::cout << ", delta^-1 = " << 1 / delta << std::endl << std::endl;

  /*
   * Вычисление наихудшего дифференцирующего вектора весов
   */
  std::cout << "Computing the worst differentiating weight vector" << std::endl;
  MatrixDn W1m = ((1 / delta) * MatrixDn::Ones(num_crit_mat, num_crit_mat)) + Calmoststar;
  auto w1 = cleany(W1m);
  auto liw1 = linearly_independent_system(w1);
  std::cout << "Linearly independent w1 = \n" << liw1.format(MatlabFmt) << std::endl << std::endl;


  /* Попробуем получить "самый наихудший" вектор.
   * Если на предыдущем шаге получился единственный вектор, то
   * данная процедура просто его нормирует.
   *
   * Иначе возьмём линейно независимые столбцы матрицы w1, отнормируем по максимуму
   * и сложим.
   */

  liw1 = liw1 / liw1.maxCoeff();
  MatrixDn worst_vector = liw1.rowwise().sum();
  std::cout << "Worst differentiating weight vector = \n" << worst_vector.format(MatlabFmt) << std::endl
            << std::endl;

  /*
   * Получим взвешенную сумму для наихудшего вектора.
   */

  MatrixDn D1 = MatrixDn::Zero(num_comp_mat, num_comp_mat);
  for (int i = 0; i < worst_vector.rows(); ++i){
      D1 += worst_vector(i) * As[i];
  }

  auto nu_1 = spectral_radius(D1);
  std::cout << "D1 = " << D1.format(MatlabFmt) << std::endl;
  std::cout << "nu_1 = " << nu_1 << std::endl;
  MatrixDn almost_x1 = (1 / nu_1) * D1;

  /*
   *  Вектор рейтингов альтернатив для D1
   */

  std::cout << "\nWorst differentiating vector of alternatives" << std::endl;

  auto x1 = cleany(almost_x1);
  auto LI_x1 = linearly_independent_system(x1);

  if (LI_x1.cols() == 1){
      // Решение единственно.

      std::cout << "x1 = \n" << LI_x1.format(MatlabFmt) << std::endl;
  } else {
      // Случай неединственности решения. Наихудший дифференцирующий вектор.

      auto delta_1 = (MatrixDn::Ones(1, x1.rows()) * x1 * MatrixDn::Ones(x1.cols(), 1)).value();
      MatrixDn almost_x1new = ((1/delta_1) * MatrixDn::Ones(x1.cols(), x1.rows())) + almost_x1;
      auto x1new = cleany(almost_x1new);

      auto LI_x1new = linearly_independent_system(x1new);

      /*
       *  Далее найдём максимальный элемент матрицы, нормируем на него,
       *  и сложим все векторы. Это и будет наихудщим вектором.
       */
      LI_x1new = LI_x1new / LI_x1new.maxCoeff();
      MatrixDn worst_x1 = LI_x1new.rowwise().sum();
      std::cout << "x1 = \n" << worst_x1.format(MatlabFmt) << std::endl;

  }

  std::cout << std::endl;

  /*
   * Вычисление наилучшего дифференцирующего вектора весов.
   * Матрица P --- линейно-независимые столбцы матрицы w.
   */
  std::cout << "\nComputing the best differenting weight vector" << std::endl;
  MatrixDn P = linearly_independent_system(w);
  std::cout << "P = \n" << P.format(MatlabFmt) << std::endl << std::endl;

  /*
   * Получение коэфициентов k и l.
   */
  auto bdfc = best_diff_vector_coefficients(P);
  std::cout << "Computing w_2" << std::endl;
  MatrixDn W2;
  for (auto& [k, l]: bdfc){
      /*
       * Вычисление w_2. Здесь получим матрицу для каждой пары (k, l),
       * объединим их в одну и найдём ЛН часть.
       */

      MatrixDn Plk_inverse = MatrixDn::Zero(P.rows(), P.cols());
      Plk_inverse(l, k) = 1 / P(l, k);
      auto t = P * (MatrixDn::Identity(P.cols(), P.cols()) +
                     Plk_inverse.transpose() * P);
      W2.conservativeResize(t.rows(), W2.cols() + t.cols());
      W2.rightCols(t.cols()) = t;
  }
  auto LI_W2 = linearly_independent_system(W2);
  std::cout << "\nLinearly independent w2 = \n" << LI_W2.format(MatlabFmt) << std::endl;

  auto LI_W2_normed = LI_W2 / LI_W2.maxCoeff();

  int LI_W2_normed_sums_min_idx;
  auto LI_W2_normed_sums = LI_W2_normed.colwise().sum().minCoeff(&LI_W2_normed_sums_min_idx);

  auto w2 = LI_W2_normed.col(LI_W2_normed_sums_min_idx) / LI_W2_normed.col(LI_W2_normed_sums_min_idx).maxCoeff();
  std::cout << "w2 = " << w2.format(MatlabFmt) << std::endl;

  MatrixDn D2 = MatrixDn::Zero(num_comp_mat, num_comp_mat);
  for (int i = 0; i < w2.rows(); ++i){
      D2 += w2(i) * As[i];
  }

  /*
   *  Вычислияем вектор рейтингов альтернатив для D_2.
   */

  std::cout << "\nComputing the best differentiating vector of alternatives" << std::endl;

  auto nu_2 = spectral_radius(D2);
  std::cout << "D2 = " << D2.format(MatlabFmt) << std::endl;
  std::cout << "nu_2 = " << nu_2 << std::endl;
  MatrixDn almost_x2 = (1 / nu_2) * D2;

  auto x2 = cleany(almost_x2);
  auto LI_x2 = linearly_independent_system(x2);

  std::cout << "\nBest differentiating vector of alternatives"  << std::endl;
  if (LI_x2.cols() == 1){
  // Случай единственности
      std::cout << "x2 = " << LI_x2.format(MatlabFmt) << std::endl;
  } else {
  // Случай неединственности
  //
  //  Результирующие столбцы матрицы нормируем по максимуму
  //  Из них выбираем тот, у которого сумма элементов минимальна
  //  То же для альтернатив
  //
      auto S = LI_x2;
      auto bdfc = best_diff_vector_coefficients(S);
      MatrixDn WS;
      for (auto& [k, l]: bdfc){
          MatrixDn Slk_inverse = MatrixDn::Zero(S.rows(), S.cols());
          Slk_inverse(l, k) = 1 / S(l, k);
          auto t = S * (MatrixDn::Identity(S.cols(), S.cols()) +
                         Slk_inverse.transpose() * S);
          WS.conservativeResize(t.rows(), WS.cols() + t.cols());
          WS.rightCols(t.cols()) = t;
      }
      auto LI_WS = linearly_independent_system(WS);
      auto LI_WS_normed = LI_WS / LI_WS.maxCoeff();

      int LI_WS_normed_sums_min_idx;
      auto LI_WS_normed_sums = LI_WS_normed.colwise().sum().minCoeff(&LI_WS_normed_sums_min_idx);

      auto x2_normed = LI_WS_normed.col(LI_WS_normed_sums_min_idx) / LI_WS_normed.col(LI_WS_normed_sums_min_idx).maxCoeff();
      std::cout << "x2 = " << x2_normed.format(MatlabFmt) << std::endl;
  }

  return 0;
}
