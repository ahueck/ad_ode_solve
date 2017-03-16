/*
 * Visualize.h
 *
 *  Created on: Feb 27, 2017
 *      Author: ahueck
 */

#ifndef VISUALIZE_H
#define VISUALIZE_H

#include <ODETypes.h>
#include <Util.h>

#include <matplotlibcpp.h>

#include <sstream>
#include <string>

namespace ode {
namespace vis {

inline void plot(const t_series& t, const y_series& y, const std::string& plot_file) {
  namespace plt = matplotlibcpp;

  const auto num_y = (*std::begin(y)).size();

  for (size_t i = 0; i < num_y; ++i) {
    scalar min = std::numeric_limits<scalar>::max();
    scalar max = std::numeric_limits<scalar>::min();
    vectory_type y_i;
    std::stringstream sstream;
    sstream << "y_" << i + 1;
    for (const auto& y_t : y) {
      const auto val = y_t[i];
      y_i.push_back(val);
      if (util::less(val, min)) {
        min = val;
      }
      if (!util::less_eq(val, max)) {
        max = val;
      }
    }
    plt::subplot(num_y, 1, i + 1);
    plt::plot(t, y_i, {{"linestyle", "-"}, {"marker", "o"}, {"color", "#006DB5" /* #0072BD */}});
    const auto ampl = (max - min) * 0.5 * 0.2;
    plt::ylim(min - ampl, max + ampl);
    plt::ylabel(sstream.str());
    plt::xlabel("t");
  }
  plt::save(plot_file);
}

} /* namespace vis */
} /* namespace ode */

#endif // VISUALIZE_H
