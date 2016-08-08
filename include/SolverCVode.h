/*
 * SolverCVode.h
 *
 *  Created on: Aug 3, 2016
 *      Author: ahueck
 */
#ifndef INCLUDE_SOLVERCVODE_H_
#define INCLUDE_SOLVERCVODE_H_

#include "Solver.h"

#include "ode.h"
#include "ode_types.h"

#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, functions, and macros */
#include <sundials/sundials_dense.h> /* definitions DlsMat and DENSE_ELEM */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

#include <assert.h>
#include <iostream>
#include <functional>

namespace ode {

namespace cvode {

class SolverCVode : public Solver<SolverCVode> {

  void* cvode_mem;
  realtype* j_buffer;

  public:
  SolverCVode() : Solver(), cvode_mem(nullptr), j_buffer(nullptr) {

  }

  SolverCVode(Eq* e, Jacobian* j) : Solver(e, j), cvode_mem(nullptr), j_buffer(nullptr) {

  }

  void f(N_Vector y, N_Vector ydot) {
    Vec_s y_v(NV_DATA_S(y), NV_LENGTH_S(y));
    Vec_s ydot_v(NV_DATA_S(ydot), NV_LENGTH_S(ydot));
    eq->f(y_v, ydot_v);
  }

  void J(N_Vector y, N_Vector fy, DlsMat J) {
    Vec_s y_v(NV_DATA_S(y), NV_LENGTH_S(y));
    // FIXME improve copying
    cvode_dense2rowmat(J);
    Mat_s J_mat(j_buffer, J->N, J->M);
    jac_f->J(y_v, J_mat, 0.0, Vec_s());
    rowmat2cvode_dense(J);
  }

  void solve();

  virtual ~SolverCVode() {
    if(cvode_mem != nullptr) {
      CVodeFree(&cvode_mem);
      cvode_mem = nullptr;
    }
    if(j_buffer != nullptr) {
      delete [] j_buffer;
      j_buffer = nullptr;
    }
  }
  private:
  void cvode_dense2rowmat(DlsMat mat) {
    assert(j_buffer != nullptr);
    const size_t n = mat->N, m = mat->M;
    for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < n; ++j) {
        j_buffer[i*n + j] = DENSE_ELEM(mat, i, j);
      }
    }
  }

  void rowmat2cvode_dense(DlsMat mat) {
    assert(j_buffer != nullptr);
    const size_t n = mat->N, m = mat->M;
    for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < n; ++j) {
        DENSE_ELEM(mat, i, j) = j_buffer[i*n + j];
      }
    }
  }
};

} /* namespace cvode */

} /* namespace ode */

extern "C" {
int f_cb(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  auto& solver = *static_cast<ode::cvode::SolverCVode*>(user_data);
  solver.f(y, ydot);
  return 0;
}

int Jac_cb(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  auto& solver = *static_cast<ode::cvode::SolverCVode*>(user_data);
  solver.J(y, fy, J);
  return 0;
}
}

void ode::cvode::SolverCVode::solve() {
  size_t NEQ = 3, OEQ = 3;

  const realtype Y1 = 1.0, Y2 = 0.0, Y3 = 0.0;
  const realtype T0 = 0.0;
  const realtype T1 = 0.1; //0.4;
  const realtype TMULT = 10.0;
  const realtype NOUT = 13.0; //12.0;
  const realtype RTOL  = 1.1e-13; //RCONST(1.0e-4)   /* scalar relative tolerance            */
  const realtype ATOL1 = 1.1e-12; //RCONST(1.0e-8)   /* vector absolute tolerance components */
  const realtype ATOL2 = 1.1e-14; //RCONST(1.0e-14)
  const realtype ATOL3 = 1.1e-12; //RCONST(1.0e-6)
  /* Create serial vector of length NEQ for I.C. and abstol */

  auto y = N_VNew_Serial(NEQ);
  auto abstol = N_VNew_Serial(NEQ);
  realtype reltol = RTOL;
  Vec_s y_w(NV_DATA_S(y), {Y1, Y2, Y3});
  Vec_s abs_w(NV_DATA_S(abstol), {ATOL1, ATOL2, ATOL3});

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  this->j_buffer = new realtype[NEQ*OEQ];

  CVodeInit(cvode_mem, f_cb, T0, y);
  void* user_data = this;
  CVodeSetUserData(cvode_mem, user_data);

  CVodeSVtolerances(cvode_mem, reltol, abstol);
  CVDense(cvode_mem, NEQ);
  CVDlsSetDenseJacFn(cvode_mem, Jac_cb);

  int flag, iout;
  realtype tout, t;
  iout = 0;  tout = T1;
  while(1) {
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    std::cout << y_w[0] << ", " << y_w[1] << ", " << y_w[2] << "\n";

    if (flag == CV_SUCCESS) {
      iout++;
      tout *= TMULT;
    }

    if (iout == NOUT) break;
  }

  /* Free y and abstol vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);
}

#endif /* INCLUDE_SOLVERCVODE_H_ */
