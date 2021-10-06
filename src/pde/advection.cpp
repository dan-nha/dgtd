#include "advection.h"

Advection::Advection(const double _flux_prefactor)
  : flux_prefactor(_flux_prefactor) {}
//------------------------------------------------------------------------
arma::mat Advection::get_initial_values(
    const arma::mat &phys_node_coords) const {

  const size_t num_rows(phys_node_coords.n_rows);
  const size_t num_cols(phys_node_coords.n_cols);
  arma::mat initial_values(num_rows, num_cols, arma::fill::zeros);
  for (size_t row(0); row < phys_node_coords.n_rows; ++row) {
    for(size_t col(0); col < phys_node_coords.n_cols; ++col) {
      initial_values(row, col) = sin(phys_node_coords(row, col));
    }
  }

  return initial_values;
}
//------------------------------------------------------------------------
std::vector<std::tuple<double, double>> 
Advection::get_surface_flux_prefactors(const size_t num_elems) const {
  
  std::vector<std::tuple<double, double>> prefactors;
  for (size_t elem(0); elem<num_elems; ++elem) {
    prefactors.push_back(
        {0.5*this->flux_prefactor, 0.5*this->flux_prefactor});
  }

  return prefactors;
}
