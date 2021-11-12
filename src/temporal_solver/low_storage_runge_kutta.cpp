#include "low_storage_runge_kutta.h"
#include "../tools/custom_errors.h"
#include "butcher_coeffs.h"

namespace TD {

Low_storage_runge_kutta::Low_storage_runge_kutta(
    const size_t _order, const size_t _num_stages)
: order(_order)
, num_stages(_num_stages)
, butcher_coeffs(this->get_butcher_coeffs(_order, _num_stages)) {
}

arma::mat Low_storage_runge_kutta::evolve_in_time(
    const std::function<arma::mat(arma::mat u, double t)> &pde,
    const arma::mat &initial_values,
    const double time,
    const double dt) const {
  
  const auto &[butcher_coeff1, butcher_coeff2, butcher_coeff3] =
      this->butcher_coeffs;

  arma::mat solution(initial_values);
  arma::mat interim_result(
      initial_values.n_rows, initial_values.n_cols, arma::fill::zeros);

  for (size_t stage(0); stage < this->num_stages; ++stage) {
    const double interim_time(time + dt * butcher_coeff3[stage]);

    interim_result *= butcher_coeff1[stage];
    interim_result += dt * pde(solution, interim_time);
    solution += butcher_coeff2[stage] * interim_result;
  }

  return solution;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
Low_storage_runge_kutta::get_butcher_coeffs(
    const size_t order,
    const size_t num_stages) const {

  if (order == 4 && num_stages == 5) {
    return {
        Lsrk_45_carpenter.butcher_coeff1,
        Lsrk_45_carpenter.butcher_coeff2,
        Lsrk_45_carpenter.butcher_coeff3};
  } else {
    throw Not_implemented(
        "Low-storage Runge-Kutta scheme of order " +
        std::to_string(order) + " and stage " +
        std::to_string(num_stages) + " currently not implemented. ");
  }
}
} // namespace TD
