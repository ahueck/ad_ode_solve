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
    ode::t_series t_single{1.7};
    ode::y_series y_single{{1.5}};
    std::vector<std::vector<ode::vectory_type>> nested_v{{{1.5, 2.5}, {3.5, 4.5}}, {{5.7, 6.7}}};
    ode::scalar sc = 0.5;
    ode::t_series t_empty;
    ode::y_series y_empty;

    std::ostringstream out;
    out << std::setprecision(1);
    out << std::fixed;

    out.str("");
    print(t_empty, out);
    REQUIRE("[]\n" == out.str());

    out.str("");
    print(y_empty, out);
    REQUIRE("[]\n" == out.str());

    out.str("");
    print(t, out);
    REQUIRE("[1.0, 2.3, 3.5, 4.0]\n" == out.str());

    out.str("");
    print(y, out);
    REQUIRE("[[1.5, 2.5], [3.0, 4.0]]\n" == out.str());

    out.str("");
    print(nested_v, out);
    REQUIRE("[[[1.5, 2.5], [3.5, 4.5]], [[5.7, 6.7]]]\n" == out.str());

    out.str("");
    print(t_single, out);
    REQUIRE("[1.7]\n" == out.str());

    out.str("");
    print(y_single, out);
    REQUIRE("[[1.5]]\n" == out.str());

    out.str("");
    print(sc, out);
    REQUIRE("0.5\n" == out.str());

    out.str("");
    print(0.5, out);
    REQUIRE("0.5\n" == out.str());

    out.str("");
    print(std::move(nested_v), out);
    REQUIRE("[[[1.5, 2.5], [3.5, 4.5]], [[5.7, 6.7]]]\n" == out.str());
  }
}
