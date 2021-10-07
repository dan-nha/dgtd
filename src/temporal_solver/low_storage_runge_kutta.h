#ifndef LOW_STORAGE_RUNGE_KUTTA_H
#define LOW_STORAGE_RUNGE_KUTTA_H

#include <armadillo>
#include <cstddef>
#include <tuple>
#include <vector>

namespace TD {
class Low_storage_runge_kutta {

public:
  Low_storage_runge_kutta(){};

  arma::mat evolve_in_time(
      const std::function<arma::mat(arma::mat u, double t)> &ode,
      const arma::mat initial_values,
      const double time,
      const double dt,
      const std::vector<double> &butcher_coeff1,
      const std::vector<double> &butcher_coeff2,
      const std::vector<double> &butcher_coeff3,
      const size_t num_stages) const;

  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
  get_butcher_coeffs(const size_t order, const size_t num_stages) const;
};
} // namespace TD

#endif
