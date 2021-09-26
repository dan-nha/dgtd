#include "pde.h"
#include <cmath>

arma::mat Pde::get_spatial_scheme(
    const arma::mat &fields,
    const double time,
    const std::vector<double> &geometric_factors,
    const arma::mat &diff_matrix,
    const arma::mat &lift_matrix,
    const double upwind_param) const {

  const arma::mat volume_fields(this->get_volume_fields(
      fields, 
      geometric_factors, 
      diff_matrix, 
      time));

  const arma::mat surface_fields(this->get_surface_fields(
      fields,
      geometric_factors,
      lift_matrix,
      upwind_param,
      time));

  return volume_fields + surface_fields;
}
//-------------------------------------------------------------------------
arma::mat Pde::get_volume_fields(
    const arma::mat &fields,
    const std::vector<double> &geometric_factors,
    const arma::mat &diff_matrix,
    const double time) const {

  const arma::mat fluxes(this->get_flux(
        fields, this->get_volume_flux_prefactor()));

  arma::mat volume_fields(-diff_matrix * fluxes);
  for (size_t elem(0); elem < fluxes.n_cols; ++elem) {
    volume_fields.col(elem) *= geometric_factors[elem];
  }

  return volume_fields;
}
//-------------------------------------------------------------------------
std::vector<std::tuple<double, double>> Pde::get_field_jumps(
    const arma::mat &fields,
    const std::tuple<double, double> boundary_conditions) const {

  auto [left_bc, right_bc] = boundary_conditions;
  const size_t last_node(fields.n_rows - 1);
  const size_t num_elems(fields.n_cols);

  std::vector<std::tuple<double, double>> field_jumps;
  if (num_elems == 1) {
    field_jumps.push_back({fields.front() - left_bc, right_bc});
  } else {
    field_jumps.push_back(
        {fields.front() - left_bc, fields(last_node, 0) - fields(0, 1)});
    if (num_elems > 2) {
      for (size_t elem(1); elem < num_elems - 1; ++elem) {
        field_jumps.push_back(
            {-std::get<1>(field_jumps.back()),
             fields(last_node, elem) - fields(0, elem + 1)});
      }
    }
    field_jumps.push_back({-std::get<1>(field_jumps.back()), right_bc});
  }

  return field_jumps;
}
//-------------------------------------------------------------------------
arma::mat Pde::get_lifted_jumps(
    const std::vector<std::tuple<double, double>> &field_jumps,
    const std::vector<std::tuple<double, double>> &flux_prefactors,
    const arma::mat &lift_matrix,
    const std::vector<double> &geometric_factors,
    const double upwind_param) const {

  const size_t num_nodes(lift_matrix.n_rows);
  const size_t num_elems(geometric_factors.size());

  arma::mat lifted_field(num_nodes, num_elems);
  for (size_t elem(0); elem < num_elems; ++elem) {

    auto [left_flux_prefactor, right_flux_prefactor] =
        flux_prefactors[elem];
    if (upwind_param == 0.) {
      const double left_face_normal(-1.);
      const double right_face_normal(1.);
      left_flux_prefactor *= left_face_normal;
      right_flux_prefactor *= right_face_normal;
    } else {
      left_flux_prefactor *= -upwind_param;
      right_flux_prefactor *= -upwind_param;
    }

    lifted_field.col(elem) = left_flux_prefactor * lift_matrix.col(0) *
                                 std::get<0>(field_jumps[elem]) +
                             right_flux_prefactor * lift_matrix.col(1) *
                                 std::get<1>(field_jumps[elem]);

    lifted_field.col(elem) *= geometric_factors[elem];
  }

  return lifted_field;
};
//-------------------------------------------------------------------------

