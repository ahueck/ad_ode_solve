/*
 * SolverCVode.cpp
 *
 *  Created on: Aug 10, 2016
 *      Author: ahueck
 */

#include <SolverCVode.h>
#include <ODE.h>
#include <ODETypes.h>

#include <nvector/nvector_serial.h>  /* serial N_Vector types, functions, and macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

#include <assert.h>
#include <iostream>
#include <functional>

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
  // FIXME improve copying
  cvode_dense2rowmat(J);
  Mat_s J_mat(j_buffer, J->N, J->M);
  jac_f->J(y_v, J_mat, 0.0, Vec_s());
  rowmat2cvode_dense(J);
}

void SolverCVode::solve(const vectory_type& y0, SolverConfig& config) {
  const size_t NEQ = config.get<int>("NEQ");
  const realtype T0 = config.get<realtype>("t0");
  const realtype TN = config.get<realtype>("tend");
  const realtype TS = config.get<realtype>("ts");
  realtype reltol = config.get<realtype>("rtol");
  auto conf_atol = config.get<std::vector<double>>("atol");
  const auto& size_atol = conf_atol.size();
  auto abstol = N_VNew_Serial(size_atol);
  for (unsigned i = 0; i < size_atol; ++i) {
    NV_Ith_S(abstol, i) = conf_atol[i];
  }

  const auto size_in = y0.size();
  auto y = N_VNew_Serial(size_in);
  for (unsigned i = 0; i < size_in; ++i) {
    NV_Ith_S(y, i) = y0[i];
  }

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

  CVodeInit(cvode_mem, cv_function_cb, T0, y);
  void* user_data = this;
  CVodeSetUserData(cvode_mem, user_data);

  CVodeSVtolerances(cvode_mem, reltol, abstol);
  CVDense(cvode_mem, NEQ);
  if (jac_f != nullptr) {
    this->j_buffer = new realtype[NEQ * size_in];
    CVDlsSetDenseJacFn(cvode_mem, cv_jacobian_dense_cb);
  }

  CVodeSetStopTime(cvode_mem, TN);
  CVodeSetMaxNumSteps(cvode_mem, -1);

  int cv_flag;
  realtype t;
  /* FIXME: output at certain intervalls */
  cv_flag = CVode(cvode_mem, TN, y, &t, CV_NORMAL);
  std::cout << t << ": " << NV_Ith_S(y, 0) << ", " << NV_Ith_S(y, 1) << ", " << NV_Ith_S(y, 2) << "\n";

  /* Free y and abstol vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);
}

inline void SolverCVode::cvode_dense2rowmat(DlsMat mat) {
  assert(j_buffer != nullptr);
  const size_t n = mat->N, m = mat->M;
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      j_buffer[i * n + j] = DENSE_ELEM(mat, i, j);
    }
  }
}

inline void SolverCVode::rowmat2cvode_dense(DlsMat mat) {
  assert(j_buffer != nullptr);
  const size_t n = mat->N, m = mat->M;
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j) {
      DENSE_ELEM(mat, i, j) = j_buffer[i * n + j];
    }
  }
}

SolverCVode::~SolverCVode() {
  if (cvode_mem != nullptr) {
    CVodeFree(&cvode_mem);
    cvode_mem = nullptr;
  }
  if (j_buffer != nullptr) {
    delete[] j_buffer;
    j_buffer = nullptr;
  }
}

} /* namespace cvode */
} /* namespace ode */

