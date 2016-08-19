#include "Rober.h"
#include "SolverCVode.h"

int main() {
  rober::Rober_s f;
  rober::ad::Rober_j j;

  std::vector<realtype> y0{1.0, 0.0, 0.0};

  ode::SolverConfig cf;
  cf.put<unsigned>("NEQ", 3u);
  cf.put<realtype>("t0", 0.0);
  cf.put<realtype>("tend", 1e11);
  cf.put<realtype>("ts", 0.1);
  cf.put<realtype>("rtol", 1.1e-13);
  std::vector<realtype> atol{1.1e-12, 1.1e-14, 1.1e-12};
  cf.put("atol", atol);

  ode::cvode::SolverCVode cv(&f, &j);
  cv.solve(y0, cf);

  return 0;
}
