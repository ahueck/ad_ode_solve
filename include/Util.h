/*
 * Util.h
 *
 *  Created on: Feb 9, 2017
 *      Author: ahueck
 */

#ifndef INCLUDE_UTIL_H_
#define INCLUDE_UTIL_H_

#include <limits>

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

} /* namespace util */
} /* namespace ode */

#endif /* INCLUDE_UTIL_H_ */
