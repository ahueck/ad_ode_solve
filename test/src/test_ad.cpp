/*
 * test_ad.cpp
 *
 *  Created on: Feb 10, 2017
 *      Author: ahueck
 */

#include <catch.hpp>

#include <UtilAD.h>

#include <TestUtil.h>

TEST_CASE("AD on a simple 1D function", "[ad_1d]") {
  using ode::test::f1d;
  using ode::test::f1d_analytical;
  using ode::test::Matrix1D;

  f1d f;
  Matrix1D J;
  J(0, 0) = 0.0;
  const std::array<double, 1> x{2.0};

  REQUIRE(J(0, 0) == Approx(0.0));

  SECTION("analytical") {
    std::array<double, 1> y{0.0};
    f(x, y);
    REQUIRE(y[0] == Approx(8.0));
    REQUIRE(f1d_analytical(2.0) == Approx(12.0));
  }

  SECTION("forward mode scalar") {
    using ode::ad::diff_fm_J;
    diff_fm_J(f, x, J, 1, 1);
    REQUIRE(J(0, 0) == Approx(12.0));
  }

  SECTION("reverse mode scalar") {
    using ode::ad::diff_rm_J;
    diff_rm_J(f, x, J, 1, 1);
    REQUIRE(J(0, 0) == Approx(12.0));
  }

  SECTION("forward mode vector") {
    using ode::ad::diff_v_fm_J;
    diff_v_fm_J<f1d, decltype(x), Matrix1D, 1, 1>(f, x, J);
    REQUIRE(J(0, 0) == Approx(12.0));
  }

  SECTION("reverse mode vector") {
    using ode::ad::diff_v_rm_J;
    diff_v_rm_J<f1d, decltype(x), Matrix1D, 1, 1>(f, x, J);
    REQUIRE(J(0, 0) == Approx(12.0));
  }
}
