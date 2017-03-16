/*
 * core_ad.h
 *
 *  Created on: Aug 8, 2016
 *      Author: ahueck
 */

#ifndef ODE_H
#define ODE_H

#include "ODETypes.h"

namespace ode {

struct Eq {
  Eq() = default;
  virtual void f(const Vec_s& y, Vec_s& dydt, const scalar t = scalar(0.0)) = 0;
  virtual ~Eq() = default;
};

struct Jacobian {
  Jacobian() = default;
  virtual void J(const Vec_s& y, Mat_s& Jf, const scalar& t = scalar(0.0), const Vec_s& dydt = Vec_s()) = 0;
  virtual ~Jacobian() = default;
};

} /* namespace ode */

#endif // ODE_H
