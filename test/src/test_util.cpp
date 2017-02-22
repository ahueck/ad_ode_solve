/*
 * test_util.cpp
 *
 *  Created on: Feb 13, 2017
 *      Author: ahueck
 */

#include <catch.hpp>

#include <Util.h>
#include <ODETypes.h>

#include <TestUtil.h>

TEST_CASE("Test general utility functions", "[utility]") {
  using namespace ode::util;

  SECTION("Comparators") {
    const double a = 1.5;
    const double b = 1.75;
    const double a_a = a;
    REQUIRE(less(a, b));
    REQUIRE(less_eq(a, b));
    REQUIRE(less_eq(a, a_a));
    REQUIRE(less_eq(a_a, a));

    REQUIRE_FALSE(less(a, a_a));
    REQUIRE_FALSE(less(b, a));
    REQUIRE_FALSE(less_eq(b, a));
  }

  SECTION("Printer") {
    ode::y_series y{{1.5, 2.5}, {3.0, 4.0}};
    ode::t_series t{1.0, 2.3, 3.5, 4.0};

    std::ostringstream out_t;
    out_t << std::setprecision(1);
    out_t << std::fixed;
    print(t, out_t);
    REQUIRE("[1.0, 2.3, 3.5, 4.0]\n" == out_t.str());

    std::ostringstream out_y;
    out_y << std::setprecision(1);
    out_y << std::fixed;
    print(y, out_y);
    REQUIRE("[[1.5, 2.5], [3.0, 4.0]]\n" == out_y.str());
  }
}
