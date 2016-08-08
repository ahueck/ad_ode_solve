/*
 * UtilAD.h
 *
 *  Created on: Aug 4, 2016
 *      Author: ahueck
 */

#ifndef INCLUDE_UTIL_AD_H_
#define INCLUDE_UTIL_AD_H_

#include <codi.hpp>

namespace ode {

namespace ad {

template<typename Function, typename Dtype, typename  Vector, typename  Matrix>
void diff_rm_J(const Function& f, const Vector& x, Matrix& J, const size_t n, const size_t m) {
  using ad_type = codi::RealReverse;
  using ad_vec = std::vector<ad_type>;

  ad_vec g_x;
  ad_vec g_f;
  for (size_t in = 0; in < n; ++in) {
    g_x.emplace_back(x[in]);
  }
  g_f.reserve(m);

  ad_type::TapeType& tape = ad_type::getGlobalTape();
  tape.setActive();
  // g_x is our current point
  for (size_t in = 0; in < n; ++in) {
    tape.registerInput(g_x[in]);
  }

  f(g_x, g_f);

  for (size_t out = 0; out < m; ++out) {
    tape.registerOutput(g_f[out]);
  }

  tape.setPassive();

  // assemble full jacobian:
  for(size_t i = 0; i < m;++i) {
    g_f[i].setGradient(1.0);
    tape.evaluate();
    for(size_t j = 0; j < n;++j) {
        J(i, j) = g_x[j].getGradient();
    }
    tape.clearAdjoints();
  }
  tape.reset();
}

} /* namespace ad */

} /* namespace ode */

#endif /* INCLUDE_UTIL_AD_H_ */
