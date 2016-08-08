/*
 * types_ode.h
 *
 *  Created on: Aug 8, 2016
 *      Author: ahueck
 */

#ifndef INCLUDE_ODE_TYPES_H_
#define INCLUDE_ODE_TYPES_H_

#include <algorithm>
#include <cstddef>
#include <initializer_list>

namespace ode {

template<typename Dtype, typename T>
struct VecWrapper {
  using value_type = Dtype;

  const size_t n;
  T v;

  VecWrapper()
    :  n{}, v{} {

  }

  VecWrapper(T vec, const size_t n)
    :  n(n), v(vec) {

  }

  VecWrapper(T vec, std::initializer_list<value_type> l)
    :  n(l.size()), v(vec) {
    std::copy(l.begin(), l.end(), &v[0]);
  }

  Dtype& operator[](const size_t i) {
    return v[i];
  }

  const Dtype& operator[](const size_t i) const {
    return v[i];
  }

  Dtype& operator()(const size_t i) {
    return v[i];
  }

  const Dtype& operator()(const size_t i) const {
    return v[i];
  }
};

template<typename Dtype, typename T>
struct MatWrapper {
  using value_type = Dtype;

  const size_t n;
  const size_t m;
  T mat;

  MatWrapper(T mat, const size_t n, const size_t m)
    : n(n), m(m), mat(mat) {

  }

  Dtype& operator()(const size_t i, const size_t j) {
    return mat[i*m + j];
  }

  const Dtype& operator()(const size_t i, const size_t j) const {
    return mat[i*m + j];
  }
};

using scalar = double;
using Vec_s = VecWrapper<scalar, scalar*>;
using Mat_s = MatWrapper<scalar, scalar*>;

} /* namespace ode */

#endif /* INCLUDE_ODE_TYPES_H_ */
