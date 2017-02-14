/*
 * test_util.cpp
 *
 *  Created on: Feb 13, 2017
 *      Author: ahueck
 */

#include <catch.hpp>

#include <Util.h>

#include <TestUtil.h>

TEST_CASE("Test general utility functions", "[utility]") {
  SECTION("Comparators") {
    using namespace ode::util;
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
}
