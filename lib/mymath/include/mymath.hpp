#pragma once

#include <ginac/ginac.h>
#include <Eigen/Core>
#include <ostream>

class MaxAlgebra {
 public:
  GiNaC::ex value;
  MaxAlgebra() : value(){};
  MaxAlgebra(GiNaC::numeric val) : value(val){};
  MaxAlgebra(int val) : value(val){};
  MaxAlgebra(GiNaC::ex val) : value(val){};
  MaxAlgebra(int num, int den) : value(GiNaC::numeric(num, den)){};

  friend std::ostream& operator<<(std::ostream& out, const MaxAlgebra& val);
  MaxAlgebra& operator+=(const MaxAlgebra& rhs);
  MaxAlgebra& operator=(const MaxAlgebra& rhs);
  MaxAlgebra& operator*=(const MaxAlgebra& rhs);
  MaxAlgebra& operator/=(const MaxAlgebra& rhs);

  explicit operator GiNaC::ex() { return value; }
};

#define OVERLOAD_OPERATOR_DECL(op, ret) \
  ret operator op(const MaxAlgebra& lhs, const MaxAlgebra& rhs);

#define OVERLOAD_OPERATOR_BOOL(op)                                 \
  bool operator op(const MaxAlgebra& lhs, const MaxAlgebra& rhs) { \
    return lhs.value op rhs.value;                                 \
  }

OVERLOAD_OPERATOR_DECL(/, MaxAlgebra);
OVERLOAD_OPERATOR_DECL(*, MaxAlgebra);

MaxAlgebra operator+(const MaxAlgebra& lhs, const MaxAlgebra& rhs);

bool isfinite(const MaxAlgebra&);

namespace Eigen {
template <>
struct NumTraits<MaxAlgebra>
    : GenericNumTraits<MaxAlgebra>  // permits to get the epsilon,
                                    // dummy_precision, lowest, highest
                                    // functions
{
  typedef MaxAlgebra Real;
  typedef MaxAlgebra NonInteger;
  typedef MaxAlgebra Nested;
  static inline int digits10() { return 0; }

  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 1,
    ReadCost = 1,
    AddCost = 3,
    MulCost = 3
  };
};

template <typename BinaryOp>
struct ScalarBinaryOpTraits<MaxAlgebra, GiNaC::ex, BinaryOp> {
  typedef MaxAlgebra ReturnType;
};

template <typename BinaryOp>
struct ScalarBinaryOpTraits<GiNaC::ex, MaxAlgebra, BinaryOp> {
  typedef MaxAlgebra ReturnType;
};
template <typename BinaryOp>
struct ScalarBinaryOpTraits<MaxAlgebra, GiNaC::numeric, BinaryOp> {
  typedef MaxAlgebra ReturnType;
};

template <typename BinaryOp>
struct ScalarBinaryOpTraits<GiNaC::numeric, MaxAlgebra, BinaryOp> {
  typedef MaxAlgebra ReturnType;
};
}  // namespace Eigen

/*
namespace GiNaC {
template<>
MaxAlgebra
pow<MaxAlgebra, MaxAlgebra>(const MaxAlgebra& a, const MaxAlgebra& b)
{}

}

*/
