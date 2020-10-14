#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include <mymath.hpp>
#include <numeric>
#include <random>
#include <vector>

#define OVERLOAD_OPERATOR_DEF(op, ret)                                         \
  ret operator op(const MaxAlgebra& lhs, const MaxAlgebra& rhs)                \
  {                                                                            \
    return lhs.value op rhs.value;                                             \
  }

OVERLOAD_OPERATOR_DEF(/, MaxAlgebra);
OVERLOAD_OPERATOR_DEF(*, MaxAlgebra);

MaxAlgebra
operator+(const MaxAlgebra& lhs, const MaxAlgebra& rhs)
{
  if (lhs.value.evalf() < rhs.value.evalf())
    return rhs;

  return lhs;
}

std::ostream&
operator<<(std::ostream& out, const MaxAlgebra& val)
{
  using std::operator<<;
  out << val.value;
  return out;
}

MaxAlgebra&
MaxAlgebra::operator=(const MaxAlgebra& rhs)
{
  this->value = rhs.value;
  return *this;
}

MaxAlgebra&
MaxAlgebra::operator+=(const MaxAlgebra& rhs)
{
  if (this->value.evalf() < rhs.value.evalf())
    this->value = rhs.value;

  return *this;
}

MaxAlgebra&
MaxAlgebra::operator*=(const MaxAlgebra& rhs)
{
  this->value *= rhs.value;
  return *this;
}

MaxAlgebra&
MaxAlgebra::operator/=(const MaxAlgebra& rhs)
{
  this->value /= rhs.value;
  return *this;
}
