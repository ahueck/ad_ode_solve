/*
 * Solver.h
 *
 *  Created on: Aug 3, 2016
 *      Author: ahueck
 */

#ifndef INCLUDE_SOLVER_H_
#define INCLUDE_SOLVER_H_

#include "ODETypes.h"

#include <map>
#include <memory>
#include <vector>

namespace ode {

class SolverConfig final {
  using value = std::shared_ptr<void>;
  using property_map = std::map<std::string, value>;

 private:
  property_map properties;

 public:
  void clear() {
    properties.clear();
  }

  bool put(const std::string& key, const char* v) {
    properties[key] = std::make_shared<std::string>(v);
    return true;
  }

  template <typename T>
  bool put(const std::string& key, const T& v) {
    properties[key] = std::make_shared<T>(v);
    return true;
  }

  template <typename T>
  const T& get(const std::string& key) const {
    auto iter_elem = properties.find(key);
    const auto val = static_cast<const T*>(iter_elem->second.get());
    return *val;
  }

  bool has(const std::string& key) const {
    return properties.find(key) != std::end(properties);
  }

  ~SolverConfig() = default;
};

struct Eq;
struct Jacobian;

template <typename Derived>
class Solver {
 protected:
  SolverConfig config;
  Eq* eq;
  Jacobian* jac_f;

 public:
  Solver() : eq(nullptr), jac_f(nullptr) {
  }

  explicit Solver(Eq* eq, Jacobian* jac_f = nullptr) : eq(eq), jac_f(jac_f) {
  }

  Derived& cast() const {
    return *static_cast<Derived*>(this);
  }

  void setConfig(const SolverConfig& config) {
    this->config = config;
  }

  void setEq(Eq* e) {
    this->eq = e;
  }

  void setJ(Jacobian* j) {
    this->jac_f = j;
  }

  std::tuple<y_series, t_series> solve(const vectory_type& y0, const SolverConfig& config) {
    return cast().solve(y0, config);
  }

  virtual ~Solver() {
  }
};

} /* namespace ode */

#endif /* INCLUDE_SOLVER_H_ */
