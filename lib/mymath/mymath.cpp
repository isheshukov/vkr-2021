#include <mymath.hpp>

#define OVERLOAD_OPERATOR_DEF(op, ret)                            \
  ret operator op(const MaxAlgebra& lhs, const MaxAlgebra& rhs) { \
    return lhs.value op rhs.value;                                \
  }

OVERLOAD_OPERATOR_DEF(/, MaxAlgebra);
OVERLOAD_OPERATOR_DEF(*, MaxAlgebra);

MaxAlgebra& abs(const MaxAlgebra& rhs) {
  auto r = GiNaC::abs(rhs.value);
  MaxAlgebra rr(r);

  return rr;
}

MaxAlgebra operator+(const MaxAlgebra& lhs, const MaxAlgebra& rhs) {
  return (lhs < rhs) ? rhs : lhs;
}

bool operator<(const MaxAlgebra& lhs, const MaxAlgebra& rhs) {
  return (lhs.value.evalf() < rhs.value.evalf()) ? true : false;
}

bool operator>(const MaxAlgebra& lhs, const MaxAlgebra& rhs) {
  return !(lhs < rhs);
}

MaxAlgebra pow(const MaxAlgebra& lhs, const MaxAlgebra& rhs) {
  return GiNaC::pow(lhs.value, rhs.value);
}

std::ostream& operator<<(std::ostream& out, const MaxAlgebra& val) {
  using std::operator<<;
  out << val.value;
  return out;
  // auto ctx = GiNaC::print_dflt(out);
  // std::stringstream buffer;
  // auto ctx = GiNaC::print_latex(buffer);
  // val.value.print(ctx);
  // auto bufStr = buffer.str();

  // bufStr.erase(remove_if(bufStr.begin(), bufStr.end(), isspace),
  // bufStr.end());

  // return out << bufStr;
}

MaxAlgebra& MaxAlgebra::operator=(const MaxAlgebra& rhs) {
  this->value = rhs.value;
  return *this;
}

MaxAlgebra& MaxAlgebra::operator+=(const MaxAlgebra& rhs) {
  if (this->value.evalf() < rhs.value.evalf()) this->value = rhs.value;

  return *this;
}

MaxAlgebra& MaxAlgebra::operator*=(const MaxAlgebra& rhs) {
  this->value *= rhs.value;
  return *this;
}

MaxAlgebra& MaxAlgebra::operator/=(const MaxAlgebra& rhs) {
  this->value /= rhs.value;
  return *this;
}
