#include "spatial_solver/elementwise_operations.h"
#include "temporal_solver/low_storage_runge_kutta.h"
#include "tools/custom_errors.h"
#include "tools/output.h"
#include "tools/get.h"

#include <typeinfo>
#include <iomanip>

namespace DGTD {
using namespace DG;
using namespace TD;

template <class Pde, class Basis, class TD_solver>
Dgtd_solver<Pde, Basis, TD_solver>::Dgtd_solver(
    Mesh::Process_mesh_data &_processed_mesh,
    const Input &_input)
    : processed_mesh(_processed_mesh),
      geo(_processed_mesh),
      input(_input),
      polynomial_order(_input.polynomial_order),
      quad_nodes(basis.get_quad_nodes(_input.polynomial_order)),
      end_time(_input.end_time), 
      dt_factor(_input.dt_factor),
      time_step(this->get_time_step()),
      diff_matrix(Elementwise_operations<Basis>(_input.polynomial_order)
                      .get_diff_matrix()),
      lift_matrix(Elementwise_operations<Basis>(_input.polynomial_order)
                      .get_lift_matrix()) {

  if (!std::is_same<TD_solver, Low_storage_runge_kutta>::value) {
    throw Not_implemented("Given time-domain solver unknown.");
  }
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
arma::mat Dgtd_solver<Pde, Basis, TD_solver>::get_solution(Pde pde) {

  const arma::mat phys_node_coords(this->get_phys_node_coords());

  Output out;
  out.store_coords(phys_node_coords);

  arma::mat fields(pde.get_initial_values(phys_node_coords));

  TD_solver lsrk(
      this->input.runge_kutta_order, this->input.runge_kutta_stages);
  
  const std::vector<double> geo_factors(this->get_geometric_factors());

  arma::mat interim_res(fields.n_rows, fields.n_cols, arma::fill::zeros);
  for (double time(0.); time <= end_time; time += this->time_step) {
    out.store_time(time);
    out.store_fields("advection", fields);
    auto dg_scheme =
        [pde, this, geo_factors] (arma::mat u, double t) {
          return pde.get_spatial_scheme(
              u,
              t,
              geo_factors,
              this->diff_matrix,
              this->lift_matrix);
        };
    fields = lsrk.evolve_in_time(
        dg_scheme,
        fields,
        time,
        this->time_step);
  }
  return fields;
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
std::vector<double> 
Dgtd_solver<Pde, Basis, TD_solver>::get_geometric_factors() {
  
  std::vector<size_t> elems;
  for (const auto region : processed_mesh.import_gmsh_regions()) {
    elems = processed_mesh.get_ordered_elems(
        processed_mesh.get_finite_elems(region));
  }

  std::vector<double> geo_factors;
  for (const auto elem: elems) {
    geo_factors.push_back(geo.get_geometric_factor(elem));
  }

  return geo_factors;
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
arma::mat
Dgtd_solver<Pde, Basis, TD_solver>::get_phys_node_coords() {

  std::vector<size_t> elems;
  for (const auto region : processed_mesh.import_gmsh_regions()) {
    elems = processed_mesh.get_ordered_elems(
        processed_mesh.get_finite_elems(region));
  }

  arma::mat phys_node_coords(this->quad_nodes.size(), elems.size());
  for (size_t n(0); n < elems.size(); ++n) {
    for (size_t m(0); m < quad_nodes.size(); ++m) {
      phys_node_coords(m, n) =
          geo.convert_ref_to_phys_coord(this->quad_nodes[m], elems[n]);
    }
  }

  return phys_node_coords;
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
double Dgtd_solver<Pde, Basis, TD_solver>::get_time_step() {

  const double min_node_dist(this->get_min_node_dist());
  const double dt(this->dt_factor*min_node_dist/(2*M_PI));
  const double num_time_steps(std::ceil(this->end_time/dt)); 
  
  return this->end_time / num_time_steps;
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
double Dgtd_solver<Pde, Basis, TD_solver>::get_min_node_dist() {

  return geo.get_min_node_dist(this->quad_nodes);
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
bool Dgtd_solver<Pde, Basis, TD_solver>::is_field_name_valid(
    const std::string &field_name) const {

  const auto field_names(Pde::get_field_names());
  if (field_names.map(field_name)) {
    return true;
  } else {
    return false;
  }
}


} // namespace DGTD
