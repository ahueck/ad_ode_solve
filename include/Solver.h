/*
 * Solver.h
 *
 *  Created on: Aug 3, 2016
 *      Author: ahueck
 */

#ifndef INCLUDE_SOLVER_H_
#define INCLUDE_SOLVER_H_

namespace ode {

struct Eq;
struct Jacobian;

template<typename Derived>
class Solver {

protected:
  Eq* eq;
  Jacobian* jac_f;

public:
  Solver() : eq(nullptr), jac_f(nullptr) {

  }

  Solver(Eq* eq, Jacobian* jac_f) : eq(eq), jac_f(jac_f) {

  }

  Derived& cast() const {
    return *static_cast<Derived*>(this);
  }

  void solve() {
    cast().solve();
  }

  virtual ~Solver() {

  }
};

} /* namespace ode */

#endif /* INCLUDE_SOLVER_H_ */
