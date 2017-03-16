/*
 * SolverCVode.cpp
 *
 *  Created on: Aug 10, 2016
 *      Author: ahueck
 */

#include <SolverCVode.h>
#include <ODE.h>
#include <ODETypes.h>
#include <Stepper.h>
#include <Util.h>
#include <UtilCVode.h>

#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, functions, and macros */
#include <sundials/sundials_math.h>  /* definition of ABS */
#include <sundials/sundials_types.h> /* definition of type realtype */

#include <cassert>
#include <functional>
#include <iostream>
#include <sstream>

namespace ode {
namespace cvode {

SolverCVode::SolverCVode() : Solver(), cvode_mem(nullptr), j_buffer(nullptr) {
}

SolverCVode::SolverCVode(Eq* e, Jacobian* j) : Solver(e, j), cvode_mem(nullptr), j_buffer(nullptr) {
}

inline void SolverCVode::f(N_Vector y, N_Vector ydot) {
  Vec_s y_v(NV_DATA_S(y), NV_LENGTH_S(y));
  Vec_s ydot_v(NV_DATA_S(ydot), NV_LENGTH_S(ydot));
  eq->f(y_v, ydot_v);
}

inline void SolverCVode::J(N_Vector y, N_Vector fy, DlsMat J) {
  Vec_s y_v(NV_DATA_S(y), NV_LENGTH_S(y));
  dlsmat2rowmat(J);
  jac_f->J(y_v, J_mat);
  rowmat2dlsmat(J);
}

std::tuple<y_series, t_series> SolverCVode::solve(const vectory_type& y0, const SolverConfig& config) {
  const bool is_stiff = config.get<bool>("stiff");
  const size_t NEQ = config.get<unsigned>("NEQ");
  const realtype T0 = config.get<realtype>("t0");
  const realtype TN = config.get<realtype>("tend");
  const realtype TS = config.get<realtype>("ts");
  realtype reltol = config.get<realtype>("rtol");

  const auto int_steps = size_t((TN - T0) / TS);  // FIXME sign etc.
  y_series y_v;
  y_v.reserve(int_steps);
  t_series t_v;
  t_v.reserve(int_steps);

  const int lmm = is_stiff ? CV_BDF : CV_ADAMS;
  const int iter = is_stiff ? CV_NEWTON : CV_FUNCTIONAL;
  auto y = container2nvector(y0);

  cvode_mem = std::unique_ptr<void, detail::cvode_del>(CVodeCreate(lmm, iter));
  auto cvode_mem_ptr = cvode_mem.get();

  CVodeInit(cvode_mem_ptr, cv_function_cb, T0, y);

  if (config.has("atolv")) {
    auto conf_atol = config.get<std::vector<realtype>>("atolv");
    auto abstol = container2nvector(conf_atol);

    CVodeSVtolerances(cvode_mem_ptr, reltol, abstol);
    N_VDestroy_Serial(abstol);
  } else {
    auto atol = config.get<realtype>("atol");
    CVodeSStolerances(cvode_mem_ptr, reltol, atol);
  }

  void* user_data = this;
  CVodeSetUserData(cvode_mem_ptr, user_data);

  CVDense(cvode_mem_ptr, NEQ);

  if (jac_f != nullptr) {
    this->j_buffer = std::make_unique<realtype[]>(NEQ * NV_LENGTH_S(y));
    J_mat.set_data(j_buffer.get());
    J_mat.set_dim(NEQ, NV_LENGTH_S(y));
    CVDlsSetDenseJacFn(cvode_mem_ptr, cv_jacobian_dense_cb);
  }

  CVodeSetStopTime(cvode_mem_ptr, TN);
  CVodeSetMaxNumSteps(cvode_mem_ptr, -1);

  auto obs = [&y_v, &t_v](N_Vector y, const realtype t) {
    y_v.push_back(nvector2container<vectory_type>(y));
    t_v.push_back(t);
  };
  auto odes = [&cvode_mem_ptr](N_Vector y, const realtype, const realtype, realtype cur_t) {
    realtype t;
    CVode(cvode_mem_ptr, cur_t, y, &t, CV_NORMAL);
  };
  ode::step_times(odes, y, T0, TN, TS, obs);

  N_VDestroy_Serial(y);

  return {y_v, t_v};
}

inline void SolverCVode::dlsmat2rowmat(DlsMat mat) {
  // N == num cols, M == num rows
  const size_t n = mat->N, m = mat->M;
  auto j_buf_ptr = j_buffer.get();
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      j_buf_ptr[i * n + j] = DENSE_ELEM(mat, i, j);
    }
  }
}

inline void SolverCVode::rowmat2dlsmat(DlsMat mat) {
  // N == num cols, M == num rows
  const size_t n = mat->N, m = mat->M;
  auto j_buf_ptr = j_buffer.get();
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      DENSE_ELEM(mat, i, j) = j_buf_ptr[i * n + j];
    }
  }
}

} /* namespace cvode */
} /* namespace ode */
