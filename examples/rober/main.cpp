#include "Rober.h"
#include "SolverCVode.h"

int main() {
  rober::Rober_s f;
  rober::ad::Rober_j j;
  std::vector<double> y0{1.0, 0.0, 0.0};
  std::vector<double> atol{1.1e-12, 1.1e-14, 1.1e-12};
  ode::SolverConfig cf;
  cf.put("NEQ", 3);
  cf.put("t0", 0.0);
  cf.put("tend", 1e11);
  cf.put("ts", 0.1);
  cf.put("rtol", 1.1e-13);
  cf.put("atol", atol);
  ode::cvode::SolverCVode cv(&f, &j);
  cv.solve(y0, cf);

  return 0;
}
