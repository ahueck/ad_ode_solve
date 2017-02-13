/*
 * TestUtil.h
 *
 *  Created on: Feb 10, 2017
 *      Author: ahueck
 */

#ifndef TEST_INCLUDE_TESTUTIL_H_
#define TEST_INCLUDE_TESTUTIL_H_

namespace ode {
namespace test {

template <size_t n>
struct Vector {
  std::array<double, n> vec;

  explicit Vector(std::initializer_list<double> v) {
    std::copy(std::begin(v), std::end(v), std::begin(vec));
  }

  double& operator()(size_t i) {
    return vec[i];
  }

  double operator()(size_t i) const {
    return vec[i];
  }

  double& operator[](size_t i) {
    return vec[i];
  }

  double operator[](size_t i) const {
    return vec[i];
  }
};

template <size_t n, size_t m>
struct Matrix {
  std::array<double, n * m> mat;

  double& operator()(size_t i, size_t j) {
    return mat[i * m + j];
  }

  double operator()(size_t i, size_t j) const {
    return mat[i * m + j];
  }

  double& operator()(size_t i) {
    return mat[i];
  }

  double operator()(size_t i) const {
    return mat[i];
  }
};

using Vector1D = Vector<1>;
using Matrix1D = Matrix<1, 1>;

} /* namespace test */
} /* namespace ode */

#endif /* TEST_INCLUDE_TESTUTIL_H_ */
