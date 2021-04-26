#pragma once

#include <ginac/ginac.h>
#include <Eigen/Core>
#include <ostream>

template<typename Op>
class MaxAlgebra {
 public:
  GiNaC::ex value;
  MaxAlgebra() : value(){};
  MaxAlgebra(GiNaC::numeric val) : value(val){};
  MaxAlgebra(int val) : value(val){};
  MaxAlgebra(GiNaC::ex val) : value(val){};
  MaxAlgebra(int num, int den) : value(GiNaC::numeric(num, den)){};

  friend MaxAlgebra<Op> operator+(const MaxAlgebra<Op>& lhs, const MaxAlgebra<Op>& rhs){
      return (lhs < rhs) ? rhs : lhs;
  };

  friend MaxAlgebra<Op> operator*(const MaxAlgebra<Op>& lhs, const MaxAlgebra<Op>& rhs);
  friend MaxAlgebra<Op> operator/(const MaxAlgebra<Op>& lhs, const MaxAlgebra<Op>& rhs);


  friend std::ostream& operator<<(std::ostream& out, const MaxAlgebra& val) {
      using std::operator<<;
      out << val.value;
      return out;
  }
  MaxAlgebra<Op>& operator=(const MaxAlgebra& rhs) {
      this->value = rhs.value;
      return *this;
  };

  MaxAlgebra<Op>& operator+=(const MaxAlgebra& rhs) {
      *this = *this + rhs;
      return *this;
  };
  MaxAlgebra<Op>& operator*=(const MaxAlgebra& rhs) {
      *this = this * rhs;
      return *this;
  };
  MaxAlgebra<Op>& operator/=(const MaxAlgebra& rhs) {
      *this = this / rhs;
      return *this;
  };

  friend bool operator<(const MaxAlgebra<Op>& lhs, const MaxAlgebra<Op>& rhs) {
      return lhs.value.evalf() < rhs.value.evalf();
  };
  friend bool operator==(const MaxAlgebra<Op>& lhs, const MaxAlgebra<Op>& rhs) {
      return lhs.value == rhs.value;
  };
  friend bool operator>(const MaxAlgebra<Op>& lhs, const MaxAlgebra<Op>& rhs) {
      return !(lhs < rhs);
  };


  MaxAlgebra<Op>& abs(const MaxAlgebra<Op>& rhs) { return MaxAlgebra<Op>(GiNaC::abs(rhs.value)); }

  explicit operator GiNaC::ex() { return value; }
  friend MaxAlgebra<Op> pow(const MaxAlgebra<Op>& lhs, const MaxAlgebra<Op>& rhs) {
      return GiNaC::pow(lhs.value, rhs.value);
  };
  friend bool isfinite(const MaxAlgebra<Op>&);
};

using MaxTimes = MaxAlgebra<std::multiplies<void>>;
using MaxPlus = MaxAlgebra<std::plus<void>>;

MaxPlus operator*(const MaxPlus& lhs, const MaxPlus& rhs) {
    return MaxPlus(lhs.value + rhs.value);
}
MaxPlus operator/(const MaxPlus& lhs, const MaxPlus& rhs) {
    return MaxPlus(lhs.value - rhs.value);
}

MaxTimes operator*(const MaxTimes& lhs, const MaxTimes& rhs) {
    return MaxTimes(lhs.value * rhs.value);
}
MaxTimes operator/(const MaxTimes& lhs, const MaxTimes& rhs) {
    return (rhs.value == 0) ? 0 : MaxTimes(lhs.value / rhs.value);
}




namespace Eigen {
template <>
struct NumTraits<MaxTimes>
    : GenericNumTraits<MaxTimes>
{
  typedef MaxTimes Real;
  typedef MaxTimes NonInteger;
  typedef MaxTimes Nested;
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
struct ScalarBinaryOpTraits<MaxTimes, GiNaC::ex, BinaryOp> {
  typedef MaxTimes ReturnType;
};

template <typename BinaryOp>
struct ScalarBinaryOpTraits<GiNaC::ex, MaxTimes, BinaryOp> {
  typedef MaxTimes ReturnType;
};
template <typename BinaryOp>
struct ScalarBinaryOpTraits<MaxTimes, GiNaC::numeric, BinaryOp> {
  typedef MaxTimes ReturnType;
};

template <typename BinaryOp>
struct ScalarBinaryOpTraits<GiNaC::numeric, MaxTimes, BinaryOp> {
  typedef MaxTimes ReturnType;
};

template <>
struct NumTraits<MaxPlus>
    : GenericNumTraits<MaxPlus>
{
  typedef MaxPlus Real;
  typedef MaxPlus NonInteger;
  typedef MaxPlus Nested;
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
struct ScalarBinaryOpTraits<MaxPlus, GiNaC::ex, BinaryOp> {
  typedef MaxTimes ReturnType;
};

template <typename BinaryOp>
struct ScalarBinaryOpTraits<GiNaC::ex, MaxPlus, BinaryOp> {
  typedef MaxTimes ReturnType;
};
template <typename BinaryOp>
struct ScalarBinaryOpTraits<MaxPlus, GiNaC::numeric, BinaryOp> {
  typedef MaxTimes ReturnType;
};

template <typename BinaryOp>
struct ScalarBinaryOpTraits<GiNaC::numeric, MaxPlus, BinaryOp> {
  typedef MaxTimes ReturnType;
};

}  // namespace Eigen

/*
namespace GiNaC {
template<>
MaxPlus
pow<MaxPlus, MaxPlus>(const MaxPlus& a, const MaxPlus& b)
{}

}

*/
