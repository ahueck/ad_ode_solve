/*
 * UtilAD.h
 *
 *  Created on: Aug 4, 2016
 *      Author: ahueck
 */

#ifndef INCLUDE_UTIL_AD_H_
#define INCLUDE_UTIL_AD_H_

#include <codi.hpp>

#include <vector>

namespace ode {

namespace ad {

template<typename Function, typename Vector, typename Matrix>
void diff_rm_J(const Function& f, const Vector& x, Matrix& J, const size_t n, const size_t m) {
  using ad_type = codi::RealReverse;
  using ad_vec = std::vector<ad_type>;

  ad_vec g_x;
  ad_vec g_f;
  for (size_t in = 0; in < n; ++in) {
    g_x.emplace_back(x[in]);
  }
  g_f.reserve(m);

  auto& tape = ad_type::getGlobalTape();
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

template<typename Function, typename Vector, typename Matrix>
void diff_fm_J(const Function& f, const Vector& x, Matrix& J, const size_t n, const size_t m) {
  using ad_type = codi::RealForward;
  using ad_vec = std::vector<ad_type>;

  ad_vec g_x;
  ad_vec g_f;
  for (size_t in = 0; in < n; ++in) {
    g_x.emplace_back(x[in]);
  }
  g_f.reserve(m);

  for (size_t i = 0; i < n; ++i) {
	g_x[i].setGradient(1.0);
	f(g_x, g_f);
	for(size_t j = 0; j < m; ++j) {
	  J(j, i) = g_f[j].getGradient();
	}
	g_x[i].setGradient(0.0);
  }
}

template<typename Function, typename Vector, typename Matrix, size_t n, size_t m>
void diff_v_fm_J(const Function& f, const Vector& x, Matrix& J) {
  using ad_type = codi::RealForwardVec<n>;

  ad_type g_x[n];
  ad_type g_f[m];
  for (size_t in = 0; in < n; ++in) {
    g_x[in] = x[in];
  }

  // FIXME pass seeding matrix
  for (size_t in = 0; in < n; ++in) {
	g_x[in].gradient()[in] = 1.0;
  }

  f(g_x, g_f);

  for(size_t j = 0; j < m; ++j) {
	for (size_t i = 0; i < n; ++i) {
	  J(j, i) = g_f[j].getGradient()[i];
	}
  }
}

template<typename Function, typename Vector, typename Matrix, size_t n, size_t m>
void diff_v_rm_J(const Function& f, const Vector& x, Matrix& J) {
  using ad_type = codi::RealReverseVec<m>;

  ad_type g_x[n];
  ad_type g_f[m];
  for (size_t in = 0; in < n; ++in) {
    g_x[in] = x[in];
  }

  auto& tape = ad_type::getGlobalTape();
  tape.setActive();

  for (size_t in = 0; in < n; ++in) {
	tape.registerInput(g_x[in]);
  }

  f(g_x, g_f);

  for (size_t out = 0; out < m; ++out) {
    tape.registerOutput(g_f[out]);
  }

  tape.setPassive();

  // FIXME pass seeding matrix
  for (size_t out = 0; out < m; ++out) {
	g_f[out].gradient()[out] = 1.0;
  }

  tape.evaluate();

  // assemble full jacobian:
  for(size_t i = 0; i < m;++i) {
    for(size_t j = 0; j < n;++j) {
        J(i, j) = g_x[j].getGradient()[i];
    }
  }
  tape.reset();
}

} /* namespace ad */

} /* namespace ode */

#endif /* INCLUDE_UTIL_AD_H_ */
