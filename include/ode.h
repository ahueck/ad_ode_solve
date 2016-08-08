/*
 * core_ad.h
 *
 *  Created on: Aug 8, 2016
 *      Author: ahueck
 */

#ifndef INCLUDE_ODE_H_
#define INCLUDE_ODE_H_

#include "ode_types.h"

namespace ode {

struct Eq {
  Eq() = default;
  virtual void f(const Vec_s& y, Vec_s& dydt, const scalar t = scalar(0.0)) = 0;
  virtual ~Eq() = default;
};

struct Jacobian {
  Jacobian() = default;
  virtual void J(const Vec_s& y, Mat_s& Jf, const scalar& t, const Vec_s& dydt) = 0;
  virtual ~Jacobian() = default;
};

} /* namespace ode */

#endif /* INCLUDE_ODE_H_ */
