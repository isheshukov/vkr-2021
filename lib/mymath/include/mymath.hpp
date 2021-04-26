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

  friend MaxAlgebra operator+(const MaxAlgebra& lhs, const MaxAlgebra& rhs){
      return (lhs < rhs) ? rhs : lhs;
  };

  friend MaxAlgebra operator*(const MaxAlgebra& lhs, const MaxAlgebra& rhs);
  friend MaxAlgebra operator/(const MaxAlgebra& lhs, const MaxAlgebra& rhs);


  friend std::ostream& operator<<(std::ostream& out, const MaxAlgebra& val) {
      using std::operator<<;
      out << val.value;
      return out;
  }
  MaxAlgebra& operator=(const MaxAlgebra& rhs) {
      this->value = rhs.value;
      return *this;
  };

  MaxAlgebra& operator+=(const MaxAlgebra& rhs) {
      *this = *this + rhs;
      return *this;
  };
  MaxAlgebra& operator*=(const MaxAlgebra& rhs) {
      *this = this * rhs;
      return *this;
  };
  MaxAlgebra& operator/=(const MaxAlgebra& rhs) {
      *this = this / rhs;
      return *this;
  };

  friend bool operator<(const MaxAlgebra& lhs, const MaxAlgebra& rhs) {
      return lhs.value.evalf() < rhs.value.evalf();
  };
  friend bool operator==(const MaxAlgebra& lhs, const MaxAlgebra& rhs) {
      return lhs.value == rhs.value;
  };
  friend bool operator>(const MaxAlgebra& lhs, const MaxAlgebra& rhs) {
      return !(lhs < rhs);
  };


  MaxAlgebra& abs(const MaxAlgebra& rhs) { return MaxAlgebra(GiNaC::abs(rhs.value)); }

  explicit operator GiNaC::ex() { return value; }
  friend MaxAlgebra pow(const MaxAlgebra& lhs, const MaxAlgebra& rhs) {
      return GiNaC::pow(lhs.value, rhs.value);
  };
  friend bool isfinite(const MaxAlgebra&);
};

MaxAlgebra<std::plus<void>> operator*(const MaxAlgebra<std::plus<void>>& lhs,
                                      const MaxAlgebra<std::plus<void>>& rhs) {
    return MaxAlgebra<std::plus<void>>(lhs.value + rhs.value);

}
MaxAlgebra<std::plus<void>> operator/(const MaxAlgebra<std::plus<void>>& lhs,
                                      const MaxAlgebra<std::plus<void>>& rhs) {
    return MaxAlgebra<std::plus<void>>(lhs.value - rhs.value);
}

MaxAlgebra<std::multiplies<void>> operator*(const MaxAlgebra<std::multiplies<void>>& lhs,
                                            const MaxAlgebra<std::multiplies<void>>& rhs) {
    return MaxAlgebra<std::multiplies<void>>(lhs.value * rhs.value);
}
MaxAlgebra<std::multiplies<void>> operator/(const MaxAlgebra<std::multiplies<void>>& lhs,
                                            const MaxAlgebra<std::multiplies<void>>& rhs) {
    return MaxAlgebra<std::multiplies<void>>(lhs.value / rhs.value);
}

using MaxTimes = MaxAlgebra<std::multiplies<void>>;
using MaxPlus = MaxAlgebra<std::plus<void>>;


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
