/*
 * test_ad.cpp
 *
 *  Created on: Feb 10, 2017
 *      Author: ahueck
 */

#include <catch.hpp>

#include <UtilAD.h>

#include <TestUtil.h>

namespace ode {
namespace test {
namespace ad {

struct f1d {
  template <typename T>
  void operator()(const T& x, T& res) const {
    res[0] = x[0] * x[0] * x[0];
  }
};

template <typename T>
inline T f1d_analytical(T x) {
  return 3 * x * x;
}

}
}
}

TEST_CASE("AD on a simple function R^1 -> R^1", "[ad_1d]") {
  using ode::test::ad::f1d;
  using ode::test::ad::f1d_analytical;
  using ode::test::Vector1D;
  using ode::test::Matrix1D;

  f1d f;
  Matrix1D J;
  J(0) = 0.0;
  const Vector1D x{2.0};

  REQUIRE(J(0, 0) == Approx(0.0));

  SECTION("analytical") {
    Vector1D y{0.0};
    f(x, y);
    REQUIRE(y[0] == Approx(8.0));
    REQUIRE(f1d_analytical(2.0) == Approx(12.0));
  }

  SECTION("forward mode scalar") {
    using ode::ad::diff_fm_J;
    diff_fm_J(f, x, J, 1, 1);
    REQUIRE(J(0) == Approx(12.0));
  }

  SECTION("reverse mode scalar") {
    using ode::ad::diff_rm_J;
    diff_rm_J(f, x, J, 1, 1);
    REQUIRE(J(0) == Approx(12.0));
  }

  SECTION("forward mode vector") {
    using ode::ad::diff_v_fm_J;
    diff_v_fm_J<f1d, Vector1D, Matrix1D, 1, 1>(f, x, J);
    REQUIRE(J(0) == Approx(12.0));
  }

  SECTION("reverse mode vector") {
    using ode::ad::diff_v_rm_J;
    diff_v_rm_J<f1d, Vector1D, Matrix1D, 1, 1>(f, x, J);
    REQUIRE(J(0) == Approx(12.0));
  }
}

