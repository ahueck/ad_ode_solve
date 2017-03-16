/*
 * Stepper.h
 *
 *  Created on: Feb 9, 2017
 *      Author: ahueck
 */

#ifndef STEPPER_H
#define STEPPER_H

#include <Util.h>

namespace ode {

template <typename ODESolver, typename State, typename Time, typename Observer>
size_t step_times(ODESolver ode_s, State state, const Time T0, const Time TE, const Time dt, Observer obs) {
  using ode::util::less_eq;
  size_t step{0};
  Time current_time_int{T0 + dt};
  obs(state, T0);
  do {
    ode_s(state, T0, TE, current_time_int);
    obs(state, current_time_int);
    ++step;
    current_time_int = T0 + (step + 1) * dt;
  } while (less_eq(current_time_int, TE));

  return step;
}

} /* namespace ode */

#endif  // STEPPER_H
