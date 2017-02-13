/*
 * test_ode.cpp
 *
 *  Created on: Feb 13, 2017
 *      Author: ahueck
 */

#include <catch.hpp>

#include <ODETypes.h>
#include <Stepper.h>

#include <TestUtil.h>

TEST_CASE("Test ode framework data type wrapper", "[ode_types]") {
  using ode::Vec_s;
  using ode::Mat_s;

  double d[2] = {1.0, 2.5};
  Vec_s v2e(d, 2);
  Mat_s m12(d, 1, 2);

  REQUIRE(v2e.n == 2);
  REQUIRE(v2e.v == d);

  REQUIRE(m12.n == 1);
  REQUIRE(m12.m == 2);
  REQUIRE(m12.mat == d);

  SECTION("Empty VectorView") {
    Vec_s v;
    REQUIRE(v.n == 0);
    REQUIRE(v.v == nullptr);
  }

  SECTION("Access VectorView") {
    REQUIRE(v2e[0] == d[0]);
    REQUIRE(v2e[1] == d[1]);
    REQUIRE(v2e(0) == d[0]);
    REQUIRE(v2e(1) == d[1]);
  }

  SECTION("Assign VectorView") {
    v2e[0] = 0.5;
    v2e[1] = 2.0;
    REQUIRE(v2e[0] == 0.5);
    REQUIRE(v2e[1] == 2.0);
    REQUIRE(d[0] == 0.5);
    REQUIRE(d[1] == 2.0);
  }

  SECTION("Empty MatrixView") {
    Mat_s m;
    REQUIRE(m.n == 0);
    REQUIRE(m.m == 0);
    REQUIRE(m.mat == nullptr);
  }

  SECTION("Access MatrixView") {
    REQUIRE(m12(0, 0) == d[0]);
    REQUIRE(m12(0, 1) == d[1]);
  }

  SECTION("Assign MatrixView") {
    m12(0, 0) = 1.5;
    m12(0, 1) = 3.0;
    REQUIRE(m12(0, 0) == 1.5);
    REQUIRE(m12(0, 1) == 3.0);
    REQUIRE(d[0] == 1.5);
    REQUIRE(d[1] == 3.0);
  }
}

TEST_CASE("Test ode framework stepper functions", "[ode_stepper]") {
  const double T0 = 0.0;
  const double TE = 100.0;
  const double dt = 10.0;
  size_t calls_solver = 0;
  size_t calls_observer = 0;
  double state = 1.0;

  auto ode_solver = [&](auto x, auto t0, auto te, auto current) {
    ++calls_solver;
    REQUIRE(t0 == Approx(T0));
    REQUIRE(te == Approx(TE));
    REQUIRE(current == Approx(T0 + dt * calls_solver));
    REQUIRE(x == Approx(1.0));
  };

  auto observer = [&](auto x, auto current) {
    ++calls_observer;
    REQUIRE(x == Approx(1.0));
  };

  auto steps = ode::step_times(ode_solver, state, T0, TE, dt, observer);

  REQUIRE(steps == calls_solver);
}
