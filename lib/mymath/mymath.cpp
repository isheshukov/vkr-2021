#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <mymath.hpp>
#include <numeric>
#include <random>
#include <vector>

#define OVERLOAD_OPERATOR(op, ret)                                             \
  ret operator op(const MaxAlgebra& lhs, const MaxAlgebra& rhs)                \
  {                                                                            \
    return MaxAlgebra(lhs.value op rhs.value);                                 \
  }

#define OVERLOAD_OPERATOR_BOOL(op)                                             \
  bool operator op(const MaxAlgebra& lhs, const MaxAlgebra& rhs)               \
  {                                                                            \
    return lhs.value op rhs.value;                                             \
  }

OVERLOAD_OPERATOR(*, MaxAlgebra)
OVERLOAD_OPERATOR(/, MaxAlgebra)

OVERLOAD_OPERATOR_BOOL(>)
OVERLOAD_OPERATOR_BOOL(<)
OVERLOAD_OPERATOR_BOOL(>=)
OVERLOAD_OPERATOR_BOOL(<=)
OVERLOAD_OPERATOR_BOOL(==)
OVERLOAD_OPERATOR_BOOL(!=)

MaxAlgebra
operator+(const MaxAlgebra& lhs, const MaxAlgebra& rhs)
{
  return MaxAlgebra(std::max(lhs.value, rhs.value));
}

std::ostream&
operator<<(std::ostream& out, const MaxAlgebra& val)
{
  using std::operator<<;
  out << val.value;
  return out;
}

MaxAlgebra&
MaxAlgebra::operator=(const MaxAlgebra rhs)
{
  this->value = rhs.value;
  return *this;
}

MaxAlgebra&
MaxAlgebra::operator+=(const MaxAlgebra rhs)
{
  this->value = (*this + rhs).value;
  return *this;
}

MaxAlgebra&
MaxAlgebra::operator*=(const MaxAlgebra rhs)
{
  value = (*this * rhs).value;
  return *this;
}

MaxAlgebra&
MaxAlgebra::operator/=(const MaxAlgebra rhs)
{
  value = (*this / rhs).value;
  return *this;
}
