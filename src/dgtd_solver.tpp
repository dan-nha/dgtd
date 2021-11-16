#include "temporal_solver/low_storage_runge_kutta.h"
#include "tools/custom_errors.h"
#include "tools/output.h"
#include "tools/get.h"

#include <typeinfo>
#include <iomanip>

namespace DGTD {
using namespace TD;

template <class Pde, class Basis, class TD_solver>
Dgtd_solver<Pde, Basis, TD_solver>::Dgtd_solver(
    Mesh::Process_mesh_data &_processed_mesh,
    const Input &_input)
    : processed_mesh{_processed_mesh},
      input{_input},
      geop(_processed_mesh),
      eop(_input.polynomial_order),
      quad_nodes{basis.get_quad_nodes(_input.polynomial_order)},
      end_time{_input.end_time}, 
      dt_factor{_input.dt_factor},
      time_step{this->get_time_step()},
      diff_matrix{eop.get_diff_matrix()},
      lift_matrix{eop.get_lift_matrix()} {

  if (!std::is_same<TD_solver, Low_storage_runge_kutta>::value) {
    throw Not_implemented("Given time-domain solver unknown.");
  }
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
arma::mat Dgtd_solver<Pde, Basis, TD_solver>::get_solution(Pde &pde) {

  const arma::mat phys_node_coords{this->get_phys_node_coords()};
  Output out;
  out.store_coords(phys_node_coords);

  TD_solver lsrk(
      this->input.runge_kutta_order, this->input.runge_kutta_stages);

  std::map<size_t, arma::mat> region_fields;
  std::map<size_t, arma::mat> region_phys_node_coords;
  std::map<size_t, std::vector<double>> region_geo_factors;
  this->initialize_dg_scheme(
      pde, region_fields, region_phys_node_coords, region_geo_factors);
  arma::mat solution(this->assemble_global_solution(
        region_fields, phys_node_coords.n_rows, phys_node_coords.n_cols));

  for (double time(0.); time <= end_time; time += this->time_step) {
    out.store_time(time);
    out.store_fields(this->input.pde_name, solution);

    for (const auto &region: this->processed_mesh.get_ordered_regions()) {
      region_fields[region] = this->evolve_dg_scheme(
            pde, 
            lsrk,
            region_fields[region], 
            time,
            region_geo_factors[region]);
    }

    solution = this->assemble_global_solution(
        region_fields, phys_node_coords.n_rows, phys_node_coords.n_cols);
  }

  return solution;
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
void Dgtd_solver<Pde, Basis, TD_solver>::initialize_dg_scheme(
    Pde &pde,
    std::map<size_t, arma::mat> &region_fields,
    std::map<size_t, arma::mat> &region_phys_node_coords,
    std::map<size_t, std::vector<double>> &region_geo_factors
    ) {
  
  for (const auto &region: this->processed_mesh.get_ordered_regions()) {
    arma::mat phys_node_coords(this->get_phys_node_coords(region));
    region_phys_node_coords[region] = phys_node_coords;
    region_fields[region] = pde.get_initial_values(phys_node_coords);
    region_geo_factors[region] = this->get_geometric_factors(region);
  }
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
arma::mat Dgtd_solver<Pde, Basis, TD_solver>::evolve_dg_scheme(
    Pde &pde,
    TD_solver &lsrk,
    const arma::mat &fields,
    const double time,
    const std::vector<double> &geo_factors) const {

    auto dg_scheme =
        [pde, this, geo_factors] (arma::mat u, double t) {
          return pde.get_spatial_scheme(
              u,
              t,
              geo_factors,
              this->diff_matrix,
              this->lift_matrix);
        };
    return lsrk.evolve_in_time(
        dg_scheme,
        fields,
        time,
        this->time_step);
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
arma::mat Dgtd_solver<Pde, Basis, TD_solver>::assemble_global_solution(
    std::map<size_t, arma::mat> &region_fields,
    const size_t num_nodes,
    const size_t num_elems) {
  arma::mat solution(num_nodes, num_elems);

  size_t global_col{0};
  for (const auto &region: this->processed_mesh.get_ordered_regions()) {
    const arma::mat region_solution(region_fields[region]);
    for (size_t col{0}; col < region_solution.n_cols; ++col, ++global_col) {
      for (size_t row{0}; row < region_solution.n_rows; ++row) {
        solution(row, global_col) = region_solution(row, col);
      }
    }
  }

  return solution;
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
std::vector<double> 
Dgtd_solver<Pde, Basis, TD_solver>::get_geometric_factors() {
  
  std::vector<size_t> elems;
  for (const auto region : processed_mesh.import_gmsh_regions()) {
    elems = this->processed_mesh.get_ordered_elems(
        this->processed_mesh.get_finite_elems(region));
  }

  std::vector<double> geo_factors;
  for (const auto elem: elems) {
    geo_factors.push_back(this->geop.get_geometric_factor(elem));
  }

  return geo_factors;
}
//---
template <class Pde, class Basis, class TD_solver>
std::vector<double> 
Dgtd_solver<Pde, Basis, TD_solver>::get_geometric_factors(
    const size_t region) {
  
  std::vector<size_t> elems;
  for (const auto region : processed_mesh.import_gmsh_regions()) {
    elems = processed_mesh.get_ordered_elems(
        processed_mesh.get_finite_elems(region));
  }

  std::vector<double> geo_factors;
  for (const auto elem: elems) {
    geo_factors.push_back(geop.get_geometric_factor(elem));
  }

  return geo_factors;
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
arma::mat
Dgtd_solver<Pde, Basis, TD_solver>::get_phys_node_coords(
    const size_t region) {

  std::vector<size_t> elems{this->processed_mesh.get_ordered_elems(
      this->processed_mesh.get_finite_elems(region))};

  arma::mat phys_node_coords(this->quad_nodes.size(), elems.size());
  for (size_t n{0}; n < elems.size(); ++n) {
    for (size_t m{0}; m < quad_nodes.size(); ++m) {
      phys_node_coords(m, n) =
          this->geop.convert_ref_to_phys_coord(this->quad_nodes[m], elems[n]);
    }
  }

  return phys_node_coords;
}
//---
template <class Pde, class Basis, class TD_solver>
arma::mat
Dgtd_solver<Pde, Basis, TD_solver>::get_phys_node_coords() {

  std::vector<size_t> elems;
  for (const auto region : processed_mesh.import_gmsh_regions()) {
    elems = processed_mesh.get_ordered_elems(
        processed_mesh.get_finite_elems(region));
  }

  arma::mat phys_node_coords(this->quad_nodes.size(), elems.size());
  for (size_t n{0}; n < elems.size(); ++n) {
    for (size_t m{0}; m < quad_nodes.size(); ++m) {
      phys_node_coords(m, n) =
          geop.convert_ref_to_phys_coord(this->quad_nodes[m], elems[n]);
    }
  }

  return phys_node_coords;
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
double Dgtd_solver<Pde, Basis, TD_solver>::get_time_step() {

  const double min_node_dist{this->get_min_node_dist()};
  const double dt{this->dt_factor*min_node_dist/(2*M_PI)};
  const double num_time_steps{std::ceil(this->end_time/dt)}; 
  
  return this->end_time / num_time_steps;
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
double Dgtd_solver<Pde, Basis, TD_solver>::get_min_node_dist() {

  return geop.get_min_node_dist(this->quad_nodes);
}
//-------------------------------------------------------------------------
template <class Pde, class Basis, class TD_solver>
bool Dgtd_solver<Pde, Basis, TD_solver>::is_field_name_valid(
    const std::string &field_name) const {

  const auto field_names{Pde::get_field_names()};
  if (field_names.map(field_name)) {
    return true;
  } else {
    return false;
  }
}

} // namespace DGTD
