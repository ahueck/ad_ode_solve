#ifndef ROBERTSON_H_
#define ROBERTSON_H_

#include <cmath>

#include "../../include/ODE.h"
#include "../../include/ODETypes.h"
#include "../../include/UtilAD.h"

namespace rober {
namespace detail {

struct rober_functor {
  template <typename Vec>
  inline void operator()(const Vec& y, Vec& dydt) const {
    using std::pow;
    dydt[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
    dydt[1] = 0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * pow(y[1], 2);
    dydt[2] = 3e7 * pow(y[1], 2);
  }
};

} /* namespace detail */

struct Rober_s : public ode::Eq {
  detail::rober_functor rober_f;

  void f(const ode::Vec_s& y, ode::Vec_s& dydt, const ode::scalar) final {
    rober_f(y, dydt);
  }
};

struct Rober_j : public ode::Jacobian {
  void J(const ode::Vec_s& y, ode::Mat_s& Jf, const ode::scalar&, const ode::Vec_s&) final {
    Jf(0, 0) = -0.04;
    Jf(0, 1) = 1e4 * y(2);
    Jf(0, 2) = 1e4 * y(1);
    Jf(1, 0) = 0.04;
    Jf(1, 1) = -1e4 * y(2) - 6e7 * y(1);
    Jf(1, 2) = -1e4 * y(1);
    Jf(2, 1) = 6e7 * y(1);
  }
};

namespace ad {

struct Rober_j : ode::Jacobian {
  detail::rober_functor rober_f;

  void J(const ode::Vec_s& y, ode::Mat_s& Jf, const ode::scalar&, const ode::Vec_s&) final {
    using namespace ode;
    ode::ad::diff_v_fm_J<detail::rober_functor, Vec_s, Mat_s, 3, 3>(rober_f, y, Jf);
  }
};

} /* namespace ad */
} /* namespace rober */

#endif
