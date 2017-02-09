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
#include <ODE.h>
#include <Stepper.h>
#include <Util.h>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <sundials/sundials_dense.h> /* definitions DlsMat and DENSE_ELEM */

#include <nvector/nvector_serial.h>  /* serial N_Vector types, functions, and macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

#include <assert.h>
#include <iostream>
#include <functional>
#include <sstream>
#include <vector>

namespace ode {

//struct MatrixView;

struct DlsMatView : public MatrixView<DlsMatView> {
  using value_type = realtype;

private:
  DlsMat matrix;

public:
  DlsMatView(const DlsMat& matrix) : matrix(matrix) {

  }

  inline value_type& operator()(const size_t i, const size_t j) {
     return DENSE_ELEM(matrix, i, j);
   }

   inline const value_type& operator()(const size_t i, const size_t j) const {
     return DENSE_ELEM(matrix, i, j);
   }
};

struct PointerCBCtr {
  std::function<void(N_Vector y, N_Vector ydot)> f;

  PointerCBCtr(std::function<void(N_Vector y, N_Vector ydot)> f) : f(f) {

  }

  template<typename... Args>
  void invoke_f(Args... args) {

  }
};

extern "C" {

inline int cv_function_cb(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
  //auto& solver = *static_cast<ode::cvode::SolverCVode*>(user_data);
  //solver.f(y, ydot);
  return 0;
}

inline int cv_jacobian_dense_cb(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void* user_data,
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  //auto& solver = *static_cast<ode::cvode::SolverCVode*>(user_data);
  //solver.J(y, fy, J);
  return 0;
}

} /* extern "C" */

namespace cvode {

template<typename Jacobian>
class SolverCVode : public Solver<SolverCVode<Jacobian>, Jacobian> {
  void* cvode_mem;
  realtype* j_buffer;

  using vectory_type = typename Solver<SolverCVode<Jacobian>, Jacobian>::vectory_type;

 public:
  SolverCVode(Eq* e, Jacobian* j) :
    Solver<SolverCVode<Jacobian>, Jacobian>(e, j),
    cvode_mem(nullptr),
    j_buffer(nullptr)
  { }

  void f(N_Vector y, N_Vector ydot) {
    Vec_s y_v(NV_DATA_S(y), NV_LENGTH_S(y));
    Vec_s ydot_v(NV_DATA_S(ydot), NV_LENGTH_S(ydot));
    //eq->f(y_v, ydot_v);
  }
  void J(N_Vector y, N_Vector fy, DlsMat J){
    Vec_s y_v(NV_DATA_S(y), NV_LENGTH_S(y));
    // FIXME improve upon copying
    cvode_dense2rowmat(J);
    Mat_s J_mat(j_buffer, J->N, J->M);
    DlsMatView view(J);
    //jac_f->J(y_v, view, 0.0, Vec_s());
    rowmat2cvode_dense(J);
  }

  void solve(const vectory_type& y0) {
    //solve(y0, config);
  }

  void solve(const vectory_type& y0, SolverConfig& config) {
    const bool is_stiff = config.get<bool>("stiff");
    const size_t NEQ = config.get<unsigned>("NEQ");
    const realtype T0 = config.get<realtype>("t0");
    const realtype TN = config.get<realtype>("tend");
    const realtype TS = config.get<realtype>("ts");
    realtype reltol = config.get<realtype>("rtol");
    auto conf_atol = config.get<std::vector<realtype>>("atol");
    auto abstol = N_VNew_Serial(conf_atol.size());
    std::copy(std::begin(conf_atol), std::end(conf_atol), NV_DATA_S(abstol));
    const auto size_in = y0.size();
    auto y = N_VNew_Serial(size_in);
    std::copy(std::begin(y0), std::end(y0), NV_DATA_S(y));

    if(is_stiff) {
      cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    } else {
      cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
    }

    CVodeInit(cvode_mem, cv_function_cb, T0, y);
    void* user_data = this;
    CVodeSetUserData(cvode_mem, user_data);

    CVodeSVtolerances(cvode_mem, reltol, abstol);
    CVDense(cvode_mem, NEQ);
    if (this->jac_f != nullptr) {
      this->j_buffer = new realtype[NEQ * size_in];
      CVDlsSetDenseJacFn(cvode_mem, cv_jacobian_dense_cb);
    }

    CVodeSetStopTime(cvode_mem, TN);
    CVodeSetMaxNumSteps(cvode_mem, -1);

    int cv_flag;
    realtype t;
    auto obs = [&NEQ](N_Vector y, const realtype t) {
      auto y_d = NV_DATA_S(y);
      std::ostringstream out_ss;
      out_ss << t << ": ";
      std::for_each(y_d, &y_d[NEQ - 1], [&] (realtype yi) {
        out_ss << yi << ", ";
      });
      out_ss << y_d[NEQ - 1] << "\n";
      std::cout << out_ss.str();
    };
    auto odes = [&] (N_Vector y, const realtype, const realtype, realtype cur_t) {
      CVode(cvode_mem, cur_t, y, &t, CV_NORMAL);
    };
    ode::step_times(odes, y, T0, TN, TS, obs);

    /* Free y and abstol vectors */
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(abstol);
  }

  virtual ~SolverCVode() {
    if (cvode_mem != nullptr) {
      CVodeFree(&cvode_mem);
      cvode_mem = nullptr;
    }
    if (j_buffer != nullptr) {
      delete[] j_buffer;
      j_buffer = nullptr;
    }
  }

 private:
  void cvode_dense2rowmat(DlsMat mat){
    // N == num cols, M == num rows
    const size_t n = mat->N, m = mat->M;
    for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < n; ++j) {
        j_buffer[i * n + j] = DENSE_ELEM(mat, i, j);
      }
    }
  }
  void rowmat2cvode_dense(DlsMat mat){
    // N == num cols, M == num rows
    const size_t n = mat->N, m = mat->M;
    for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < n; ++j) {
        DENSE_ELEM(mat, i, j) = j_buffer[i * n + j];
      }
    }
  }
};

} /* namespace cvode */
} /* namespace ode */


#endif /* INCLUDE_SOLVERCVODE_H_ */
