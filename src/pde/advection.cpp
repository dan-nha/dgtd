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
arma::mat Advection::get_surface_fields(
      const arma::mat &fields,
      const std::vector<double> &geometric_factors,
      const arma::mat &lift_matrix,
      const double upwind_param,
      const double time) const {

  std::tuple<double, double> boundary_conditions(
      this->get_boundary_conditions(fields, time));
  std::vector<std::tuple<double, double>> field_jumps(
      Pde::get_field_jumps(fields, boundary_conditions));

  const size_t num_elems(fields.n_cols);
  std::vector<std::tuple<double, double>> surface_flux_prefactors(
      this->get_surface_flux_prefactors(num_elems));

  return get_lifted_jumps(
      field_jumps,
      surface_flux_prefactors,
      lift_matrix,
      geometric_factors,
      upwind_param);
}
//------------------------------------------------------------------------
std::vector<std::tuple<double, double>> 
Advection::get_surface_flux_prefactors(const size_t num_elems) const {
  
  std::vector<std::tuple<double, double>> prefactors;
  for (size_t elem(0); elem<num_elems; ++elem) {
    prefactors.push_back({this->flux_prefactor, this->flux_prefactor});
  }

  return prefactors;
}
