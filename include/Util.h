/*
 * Util.h
 *
 *  Created on: Feb 9, 2017
 *      Author: ahueck
 */

#ifndef INCLUDE_UTIL_H_
#define INCLUDE_UTIL_H_

#include <iostream>
#include <limits>
#include <sstream>

namespace ode {
namespace util {

template <typename T>
inline bool less(T a, T b) {
  return b - a > std::numeric_limits<T>::epsilon();
}

template <typename T>
inline bool less_eq(T a, T b) {
  return a - b <= std::numeric_limits<T>::epsilon();
}

namespace detail {

template <typename T>
using is_scalar_t = typename std::is_scalar<T>::type;

template <typename T>
inline void print(std::true_type, const T& v, std::ostream& out) {
  out << v;
}

template <typename Container>
inline void print(std::false_type, const Container& v, std::ostream& out) {
  out << "[";
  auto last = std::end(v) - 1;
  std::for_each(std::begin(v), last, [&out](const auto& yi) {
    detail::print(is_scalar_t<std::remove_reference_t<decltype(yi)>>(), yi, out);
    out << ", ";
  });
  detail::print(is_scalar_t<std::remove_reference_t<decltype(*last)>>(), *last, out);
  out << "]";
}
}

template <typename T>
inline void print(const T& v, std::ostream& out = std::cout) {
  using detail::is_scalar_t;
  detail::print(is_scalar_t<T>(), v, out);
  out << std::endl;
}

} /* namespace util */
} /* namespace ode */

#endif /* INCLUDE_UTIL_H_ */
