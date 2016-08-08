#include "ode_systems/rober.h"
#include "SolverCVode.h"

int main() {
  rober::Rober_s f;
  rober::ad::Rober_j j;
  ode::cvode::SolverCVode cv(&f, &j);
  cv.solve();

  return 0;
}
