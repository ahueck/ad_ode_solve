/*
 * SolverCVode.h
 *
 *  Created on: Aug 3, 2016
 *      Author: ahueck
 */
#ifndef INCLUDE_SOLVERCVODE_H_
#define INCLUDE_SOLVERCVODE_H_

#include "Solver.h"

#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <sundials/sundials_dense.h> /* definitions DlsMat and DENSE_ELEM */

namespace ode {

namespace cvode {

class SolverCVode : public Solver<SolverCVode> {

  void* cvode_mem;
  realtype* j_buffer;

  public:
  SolverCVode();
  SolverCVode(Eq* e, Jacobian* j);

  void f(N_Vector y, N_Vector ydot);
  void J(N_Vector y, N_Vector fy, DlsMat J);
  void solve();

  virtual ~SolverCVode();

  private:

  void cvode_dense2rowmat(DlsMat mat);
  void rowmat2cvode_dense(DlsMat mat);
};

} /* namespace cvode */

} /* namespace ode */

extern "C" {

inline int cv_function_cb(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  auto& solver = *static_cast<ode::cvode::SolverCVode*>(user_data);
  solver.f(y, ydot);
  return 0;
}

inline int cv_jacobian_dense_cb(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  auto& solver = *static_cast<ode::cvode::SolverCVode*>(user_data);
  solver.J(y, fy, J);
  return 0;
}

} /* extern "C" */

#endif /* INCLUDE_SOLVERCVODE_H_ */
