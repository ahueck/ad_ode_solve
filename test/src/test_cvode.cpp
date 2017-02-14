/*
 * test_cvode.cpp
 *
 *  Created on: Feb 14, 2017
 *      Author: ahueck
 */

#include <catch.hpp>

#include <UtilCVode.h>

#include <TestUtil.h>

TEST_CASE("Test general CVode utility functions", "[cvode_utility]") {
  using namespace ode::cvode;
  std::vector<double> d{1.0, 2.0};

  SECTION("Container to nvector") {
    auto nv = container2nvector(d);
    REQUIRE(NV_Ith_S(nv, 0) == d[0]);
    REQUIRE(NV_Ith_S(nv, 1) == d[1]);
    auto data = NV_DATA_S(nv);
    REQUIRE_FALSE(data == d.data());
    data[0] = 2.5;
    REQUIRE_FALSE(NV_Ith_S(nv, 0) == d[0]);
    REQUIRE(NV_Ith_S(nv, 0) == 2.5);
  }

  SECTION("nvector to container") {
    auto nv = container2nvector(d);
    REQUIRE(NV_Ith_S(nv, 0) == d[0]);
    REQUIRE(NV_Ith_S(nv, 1) == d[1]);

    auto vec = nvector2container<std::vector<double>>(nv);
    REQUIRE(NV_Ith_S(nv, 0) == vec[0]);
    REQUIRE(NV_Ith_S(nv, 1) == vec[1]);
  }
}
