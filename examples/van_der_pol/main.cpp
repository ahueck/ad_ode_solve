/*
 * main.cpp
 *
 *  Created on: Feb 10, 2017
 *      Author: ahueck
 */

#include "VDP.h"

#include "SolverCVode.h"
#include "Visualize.h"

#include "cxxopts.hpp"

int main(int argc, char** argv) {
  using ode::y_series;
  using ode::t_series;
  namespace opt = cxxopts;

  opt::Options options("Van-der-Pol ODE", "Solution to the stiff Van der Pol ODE");
  options.add_options()("h,help", "Print help");
  options.add_options("Time")("t0", "Start time", opt::value<realtype>()->default_value("0.0"))(
      "tN", "End time", opt::value<realtype>()->default_value("3000.0"))("ts", "Time step",
                                                                         opt::value<realtype>()->default_value("10.0"));

  options.add_options("Tolerance")("r,relative", "Relative tolerances",
                                   opt::value<realtype>()->default_value("1.1e-13"))(
      "a,absolute", "Absolute tolerances", opt::value<realtype>()->default_value("1.0e-6"))(
      "absolute-vector", "Absolute tolerances vector", opt::value<std::vector<realtype>>());

  options.add_options("Initial value")("i,initial", "Initial value of ODE", opt::value<std::vector<realtype>>());

  options.parse(argc, argv);

  if (options.count("help") > 0) {
    std::cout << options.help({"", "Time", "Tolerance", "Initial value"});
  }

  vdp::VDP_s f;
  vdp::ad::VDP_j j;
  auto y0 = options.count("initial") > 0 ? options["i"].as<std::vector<realtype>>() : std::vector<realtype>{2.0, 0.0};

  ode::SolverConfig cf;
  cf.put<bool>("stiff", true);
  cf.put<unsigned>("NEQ", 2u);
  cf.put<realtype>("t0", options["t0"].as<realtype>());
  cf.put<realtype>("tend", options["tN"].as<realtype>());
  cf.put<realtype>("ts", options["ts"].as<realtype>());
  cf.put<realtype>("rtol", options["r"].as<realtype>());
  cf.put("atol", options["r"].as<realtype>());

  ode::cvode::SolverCVode cv(&f, &j);
  y_series y;
  t_series t;
  std::tie(y, t) = cv.solve(y0, cf);

  ode::vis::plot(t, y, "./vdp.png");

  return 0;
}
