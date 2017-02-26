/*
 * UtilCVode.h
 *
 *  Created on: Feb 10, 2017
 *      Author: ahueck
 */

#ifndef INCLUDE_UTILCVODE_H_
#define INCLUDE_UTILCVODE_H_

namespace ode {
namespace cvode {

#include <nvector/nvector_serial.h> /* serial N_Vector types, functions, and macros */

#include <algorithm>

template <typename Container>
inline N_Vector container2nvector(const Container& c) {
  const auto s = std::begin(c);
  const auto e = std::end(c);
  auto vec = N_VNew_Serial(std::distance(s, e));
  std::copy(s, e, NV_DATA_S(vec));
  return vec;
}

template <typename Container>
inline Container nvector2container(N_Vector c) {
  Container vec;
  std::copy(NV_DATA_S(c), NV_DATA_S(c) + NV_LENGTH_S(c), std::back_inserter(vec));
  return vec;
}

} /* namespace cvode */
} /* namespace ode */

#endif /* INCLUDE_UTILCVODE_H_ */
