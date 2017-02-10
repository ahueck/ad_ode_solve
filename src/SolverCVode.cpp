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

#include <nvector/nvector_serial.h>  /* serial N_Vector types, functions, and macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

#include <assert.h>
#include <iostream>
#include <functional>
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
  // FIXME improve upon copying
  dlsmat2rowmat(J);
  Mat_s J_mat(j_buffer, J->N, J->M);
  jac_f->J(y_v, J_mat, 0.0, Vec_s());
  rowmat2dlsmat(J);
}

SolverCVode::vectory_type SolverCVode::solve(const vectory_type& y0, SolverConfig& config) {
  const bool is_stiff = config.get<bool>("stiff");
  const size_t NEQ = config.get<unsigned>("NEQ");
  const realtype T0 = config.get<realtype>("t0");
  const realtype TN = config.get<realtype>("tend");
  const realtype TS = config.get<realtype>("ts");
  realtype reltol = config.get<realtype>("rtol");

  const int lmm = is_stiff ? CV_BDF : CV_ADAMS;
  const int iter = is_stiff ? CV_NEWTON : CV_FUNCTIONAL;
  auto y = container2nvector(y0);

  cvode_mem = std::unique_ptr<void, detail::cvode_del>(CVodeCreate(lmm, iter));
  auto cvode_mem_ptr = cvode_mem.get();

  CVodeInit(cvode_mem_ptr, cv_function_cb, T0, y);

  if(config.has("atolv")) {
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
    this->j_buffer = new realtype[NEQ * NV_LENGTH_S(y)];
    CVDlsSetDenseJacFn(cvode_mem_ptr, cv_jacobian_dense_cb);
  }

  CVodeSetStopTime(cvode_mem_ptr, TN);
  CVodeSetMaxNumSteps(cvode_mem_ptr, -1);

  int cv_flag;
  realtype t;
  auto obs = [&NEQ](N_Vector y, const realtype t) {
    auto y_d = NV_DATA_S(y);
    std::ostringstream out_ss;
    out_ss << t << ": ";
    std::for_each(y_d, &y_d[NEQ - 1], [&y_d, &out_ss] (realtype yi) {
      out_ss << yi << ", ";
    });
    out_ss << y_d[NEQ - 1] << "\n";
    std::cout << out_ss.str();
  };
  auto odes = [&] (N_Vector y, const realtype, const realtype, realtype cur_t) {
    CVode(cvode_mem_ptr, cur_t, y, &t, CV_NORMAL);
  };
  ode::step_times(odes, y, T0, TN, TS, obs);

  auto y_N = nvector2container<vectory_type>(y);
  N_VDestroy_Serial(y);

  return y_N;
}

inline void SolverCVode::dlsmat2rowmat(DlsMat mat) {
  // N == num cols, M == num rows
  const size_t n = mat->N, m = mat->M;
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      j_buffer[i * n + j] = DENSE_ELEM(mat, i, j);
    }
  }
}

inline void SolverCVode::rowmat2dlsmat(DlsMat mat) {
  // N == num cols, M == num rows
  const size_t n = mat->N, m = mat->M;
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      DENSE_ELEM(mat, i, j) = j_buffer[i * n + j];
    }
  }
}

SolverCVode::~SolverCVode() {
  if (j_buffer != nullptr) {
    delete[] j_buffer;
    j_buffer = nullptr;
  }
}

} /* namespace cvode */
} /* namespace ode */

