/*
 * VDP.h
 *
 *  Created on: Feb 10, 2017
 *      Author: ahueck
 */

#ifndef VAN_DER_POL_VDP_H_
#define VAN_DER_POL_VDP_H_


#include <cmath>

#include "../../include/ODE.h"
#include "../../include/ODETypes.h"
#include "../../include/UtilAD.h"

namespace vdp {

static const auto mu = 1000.0;

namespace detail {

struct vdp_functor {
  template <typename Vec>
  inline void operator()(const Vec& y, Vec& dydt) const {
    using std::pow;
    dydt[0] = y[1];
    dydt[1] = -y[0] - mu * y[1] * (pow(y[0], 2) - 1.0);
  }
};

} /* namespace detail */

struct VDP_s : public ode::Eq {
  detail::vdp_functor vdp_f;

  void f(const ode::Vec_s& y, ode::Vec_s& dydt, const ode::scalar) final {
    vdp_f(y, dydt);
  }
};

struct VDP_j : public ode::Jacobian {
  void J(const ode::Vec_s& y, ode::Mat_s& Jf, const ode::scalar&, const ode::Vec_s&) final {
    using std::pow;
    Jf(0, 0) = 0.0;
    Jf(0, 1) = 1.0;
    Jf(1, 0) = -1.0 - 2.0 * mu * y(0) * y(1);
    Jf(1, 1) = -mu * (pow(y(0), 2) - 1.0);
  }
};

namespace ad {

struct VDP_j : ode::Jacobian {
  detail::vdp_functor vdp_f;

  void J(const ode::Vec_s& y, ode::Mat_s& Jf, const ode::scalar&, const ode::Vec_s&) final {
    using namespace ode;
    ode::ad::diff_v_fm_J<detail::vdp_functor, Vec_s, Mat_s, 2, 2>(vdp_f, y, Jf);
  }
};

} /* namespace ad */
} /* namespace vdp */

#endif /* VAN_DER_POL_VDP_H_ */
