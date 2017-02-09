/*
 * Stepper.h
 *
 *  Created on: Feb 9, 2017
 *      Author: ahueck
 */

#ifndef INCLUDE_STEPPER_H_
#define INCLUDE_STEPPER_H_

#include <Util.h>

namespace ode {

template<typename ODESolver, typename State, typename Time, typename Observer>
size_t step_times(ODESolver odes, State state, const Time T0, const Time TE, const Time dt, Observer obs) {
  using ode::util::less;
  size_t step{0};
  Time current_time_int{T0 + dt};
  obs(state, T0);
  do {
    odes(state, T0, TE, current_time_int);
    obs(state, current_time_int);
    ++step;
    current_time_int = T0 + (step + 1) * dt;
  } while(less(current_time_int, TE));

  return step;
}

} /* namespace ode */

#endif /* INCLUDE_STEPPER_H_ */
