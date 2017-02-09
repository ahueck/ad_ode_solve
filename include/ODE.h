/*
 * core_ad.h
 *
 *  Created on: Aug 8, 2016
 *      Author: ahueck
 */

#ifndef INCLUDE_ODE_H_
#define INCLUDE_ODE_H_

#include "ODETypes.h"

namespace ode {

struct Eq {
  Eq() = default;
  virtual void f(const Vec_s& y, Vec_s& dydt, const scalar t = scalar(0.0)) = 0;
  virtual ~Eq() = default;
};

template <typename Jac>
struct Jacobian {
  Jacobian() = default;

  Jac& cast() {
    return *static_cast<Jac*>(this);
  }

  template <typename Matrix>
  void J(const Vec_s& y, MatrixView<Matrix>& Jf, const scalar& t, const Vec_s& dydt) {
    cast().J(y, Jf, t, dydt);
  }

  virtual ~Jacobian() = default;
};

} /* namespace ode */

#endif /* INCLUDE_ODE_H_ */
