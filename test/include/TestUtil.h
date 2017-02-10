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

struct f1d {
  template <typename T>
  void operator()(const T& x, T& res) const {
    res[0] = x[0] * x[0] * x[0];
  }
};

template <typename T>
inline T f1d_analytical(T x) {
  return 3 * x * x;
}

template<size_t n, size_t m>
struct Matrix {
  std::array<double, n*m> mat;

  double& operator()(size_t i, size_t j) {
    return mat[i*m + j];
  }

  double operator()(size_t i, size_t j) const {
    return mat[i*m + j];
  }
};

using Matrix1D = Matrix<1,1>;

} /* namespace test */
} /* namespace ode */

#endif /* TEST_INCLUDE_TESTUTIL_H_ */
