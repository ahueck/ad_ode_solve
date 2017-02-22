/*
 * types_ode.h
 *
 *  Created on: Aug 8, 2016
 *      Author: ahueck
 */

#ifndef INCLUDE_ODETYPES_H_
#define INCLUDE_ODETYPES_H_

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <vector>

namespace ode {

template <typename T>
struct VectorView {
  using value_type = T;
  using pointer_type = std::add_pointer_t<value_type>;

  pointer_type v;
  const size_t n;

  VectorView() : v(nullptr), n(0) {
  }

  VectorView(pointer_type vec, const size_t n) : v(vec), n(n) {
  }

  VectorView(pointer_type vec, std::initializer_list<value_type> l) : v(vec), n(l.size()) {
    std::copy(std::begin(l), std::end(l), &v[0]);
  }

  void set_data(pointer_type v) {
    this->v = v;
  }

  void set_dim(size_t n) {
    this->n = n;
  }

  inline value_type& operator[](const size_t i) {
    return v[i];
  }

  inline const value_type& operator[](const size_t i) const {
    return v[i];
  }

  inline value_type& operator()(const size_t i) {
    return v[i];
  }

  inline const value_type& operator()(const size_t i) const {
    return v[i];
  }
};

template <typename T>
struct MatrixView {
  using value_type = T;
  using pointer_type = std::add_pointer_t<value_type>;

  pointer_type mat;
  size_t n;  // #columns
  size_t m;  // #rows

  MatrixView() : mat(nullptr), n(0), m(0) {
  }

  MatrixView(pointer_type mat, const size_t n, const size_t m) : mat(mat), n(n), m(m) {
  }

  void set_data(pointer_type mat) {
    this->mat = mat;
  }

  void set_dim(size_t cols, size_t rows) {
    this->n = cols;
    this->m = rows;
  }

  inline value_type& operator()(const size_t i, const size_t j) {
    return mat[i * n + j];
  }

  inline const value_type& operator()(const size_t i, const size_t j) const {
    return mat[i * n + j];
  }
};

using scalar = double;
using Vec_s = VectorView<scalar>;
using Mat_s = MatrixView<scalar>;
using vectory_type = std::vector<scalar>;
using y_series = std::vector<vectory_type>;
using t_series = vectory_type;

} /* namespace ode */

#endif /* INCLUDE_ODETYPES_H_ */
