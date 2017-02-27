/*
 * main.cpp
 *
 *  Created on: Feb 10, 2017
 *      Author: ahueck
 */

#include "VDP.h"
#include "SolverCVode.h"
#include "Visualize.h"

int main() {
  using ode::y_series;
  using ode::t_series;

  vdp::VDP_s f;
  vdp::ad::VDP_j j;
  std::vector<realtype> y0{2.0, 0.0};

  ode::SolverConfig cf;
  cf.put<bool>("stiff", true);
  cf.put<unsigned>("NEQ", 2u);
  cf.put<realtype>("t0", 0.0);
  cf.put<realtype>("tend", 3e3);
  cf.put<realtype>("ts", 1e1);
  cf.put<realtype>("rtol", 1.1e-13);
  cf.put("atol", 1.0e-6);

  ode::cvode::SolverCVode cv(&f, &j);
  y_series y;
  t_series t;
  std::tie(y, t) = cv.solve(y0, cf);

  ode::vis::plot(t, y, "./vdp.png");

  return 0;
}
