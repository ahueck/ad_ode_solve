/*
 * SolverCVode.h
 *
 *  Created on: Aug 3, 2016
 *      Author: ahueck
 */
#ifndef INCLUDE_SOLVERCVODE_H_
#define INCLUDE_SOLVERCVODE_H_

#include "Solver.h"
#include "ODETypes.h"

#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <sundials/sundials_dense.h> /* definitions DlsMat and DENSE_ELEM */

#include <memory>

namespace ode {

namespace cvode {

namespace detail {

struct cvode_del final {
  void operator()(void* cvode_mem) const {
    CVodeFree(&cvode_mem);
  }
};

} /* namespace detail */

class SolverCVode : public Solver<SolverCVode> {
  std::unique_ptr<void, detail::cvode_del> cvode_mem;
  std::unique_ptr<realtype[]> j_buffer;
  Mat_s J_mat;

 public:
  SolverCVode();
  SolverCVode(Eq* e, Jacobian* j);

  void f(N_Vector y, N_Vector ydot);
  void J(N_Vector y, N_Vector fy, DlsMat J);

  Solver<SolverCVode>::vectory_type solve(const vectory_type& y0, SolverConfig& config);

  virtual ~SolverCVode() = default;

 private:
  void dlsmat2rowmat(DlsMat mat);
  void rowmat2dlsmat(DlsMat mat);
};

} /* namespace cvode */
} /* namespace ode */

extern "C" {

inline int cv_function_cb(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
  auto& solver = *static_cast<ode::cvode::SolverCVode*>(user_data);
  solver.f(y, ydot);
  return 0;
}

inline int cv_jacobian_dense_cb(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void* user_data,
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  auto& solver = *static_cast<ode::cvode::SolverCVode*>(user_data);
  solver.J(y, fy, J);
  return 0;
}

} /* extern "C" */

#endif /* INCLUDE_SOLVERCVODE_H_ */
