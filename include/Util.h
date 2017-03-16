/*
 * Util.h
 *
 *  Created on: Feb 9, 2017
 *      Author: ahueck
 */

#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <limits>
#include <sstream>

namespace ode {
namespace util {

namespace detail {

template <typename T>
using is_scalar_t = typename std::is_scalar<T>::type;

template <typename T>
using is_scalar_rem_ref_t = is_scalar_t<std::remove_reference_t<T>>;

template <typename T>
inline void print(std::true_type, const T& v, std::ostream& out) {
  out << v;
}

template <typename Container>
inline void print(std::false_type, Container&& v, std::ostream& out) {
  out << "[";

  const auto last = std::prev(std::end(v));
  const auto begin = std::begin(v);
  const auto dist = std::distance(begin, last);
  if (dist >= 0) {
    std::for_each(begin, last, [&out](const auto& yi) {
      detail::print(is_scalar_rem_ref_t<decltype(yi)>(), yi, out);
      out << ", ";
    });
    detail::print(is_scalar_rem_ref_t<decltype(*last)>(), *last, out);
  }

  out << "]";
}

} /* namespace detail */

template <typename T>
inline bool less(const T& a, const T& b) {
  return b - a > std::numeric_limits<T>::epsilon();
}

template <typename T>
inline bool less_eq(const T& a, const T& b) {
  return a - b <= std::numeric_limits<T>::epsilon();
}

template <typename T>
inline bool eq(const T& a, const T& b) {
  return a - b == std::numeric_limits<T>::epsilon();
}

template <typename T>
inline void print(T&& v, std::ostream& out = std::cout) {
  using detail::is_scalar_rem_ref_t;
  detail::print(is_scalar_rem_ref_t<T>(), std::forward<T>(v), out);
  out << std::endl;
}

} /* namespace util */
} /* namespace ode */

#endif  // UTIL_H
