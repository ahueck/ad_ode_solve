/*
 * Solver.h
 *
 *  Created on: Aug 3, 2016
 *      Author: ahueck
 */

#ifndef INCLUDE_SOLVER_H_
#define INCLUDE_SOLVER_H_

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
  SolverConfig() = default;

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
  const T& get(const std::string& key) {
    auto val = static_cast<T*>(properties[key].get());
    return *val;
  }
};

struct Eq;
//struct Jacobian;

template <typename Derived, typename Jacobian>
class Solver {
 public:
  SolverConfig config;
  Eq* eq;
  Jacobian* jac_f;

 public:
  using vectory_type = std::vector<double>;

  Solver(Eq* eq, Jacobian* j) : eq(eq), jac_f(j) {

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

  void solve(const vectory_type& y0) {
    cast().solve(y0);
  }

  void solve(const vectory_type& y0, SolverConfig& config) {
    cast().solve(y0, config);
  }

  virtual ~Solver() {
  }
};

/*
template <typename Derived>
struct PointerCtr {

  Derived& cast() const {
    return *static_cast<Derived*>(this);
  }

  template<typename... Args>
  void invoke_f(Args... args) {
    cast().f(args);
  }
};
*/

} /* namespace ode */

#endif /* INCLUDE_SOLVER_H_ */
