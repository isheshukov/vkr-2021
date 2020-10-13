#pragma once

#include <Eigen/Core>
#include <ginac/ginac.h>
#include <ostream>

using namespace GiNaC;

class MaxAlgebra
{
public:
  numeric value;
  MaxAlgebra()
    : value(){};
  MaxAlgebra(numeric val)
    : value(val){};
  MaxAlgebra(int num, int den)
    : value(num, den){};

  friend std::ostream& operator<<(std::ostream& out, const MaxAlgebra& val);
  MaxAlgebra& operator+=(const MaxAlgebra& rhs);
  MaxAlgebra operator=(const MaxAlgebra& rhs);
  MaxAlgebra operator*=(const MaxAlgebra& rhs);
  MaxAlgebra operator/=(const MaxAlgebra& rhs);

  explicit operator numeric() { return value; }
};

MaxAlgebra
operator+(const MaxAlgebra& lhs, const MaxAlgebra& rhs);
MaxAlgebra
operator*(const MaxAlgebra& lhs, const MaxAlgebra& rhs);

bool
isfinite(const MaxAlgebra&);

namespace Eigen {
template<>
struct NumTraits<MaxAlgebra>
  : GenericNumTraits<MaxAlgebra> // permits to get the epsilon, dummy_precision,
                                 // lowest, highest functions
{
  typedef MaxAlgebra Real;
  typedef MaxAlgebra NonInteger;
  typedef MaxAlgebra Nested;
  static inline int digits10() { return 0; }

  enum
  {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};

template<typename BinaryOp>
struct ScalarBinaryOpTraits<MaxAlgebra, numeric, BinaryOp>
{
  typedef MaxAlgebra ReturnType;
};

template<typename BinaryOp>
struct ScalarBinaryOpTraits<numeric, MaxAlgebra, BinaryOp>
{
  typedef MaxAlgebra ReturnType;
};
}
