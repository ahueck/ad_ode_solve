#include "Rober.h"
#include "SolverCVode.h"
#include "Util.h"
#include "Visualize.h"

int main() {
  using ode::y_series;
  using ode::t_series;

  rober::Rober_s f;
  rober::ad::Rober_j j;
  std::vector<realtype> y0{1.0, 0.0, 0.0};

  ode::SolverConfig cf;
  cf.put<bool>("stiff", true);
  cf.put<unsigned>("NEQ", 3u);
  cf.put<realtype>("t0", 0.0);
  cf.put<realtype>("tend", 1e11);
  cf.put<realtype>("ts", 1e9);
  cf.put<realtype>("rtol", 1.1e-13);
  cf.put<std::vector<realtype>>("atolv", {1.1e-12, 1.1e-14, 1.1e-12});

  ode::cvode::SolverCVode cv(&f, &j);
  y_series y;
  t_series t;
  std::tie(y, t) = cv.solve(y0, cf);

  ode::vis::plot(t, y, "./rober.png");

  return 0;
}
