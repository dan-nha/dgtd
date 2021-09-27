#include "spatial_solver/elementwise_operations.h"
#include "spatial_solver/geometric_operations.h"
#include "temporal_solver/low_storage_runge_kutta.h"
#include "spatial_solver/mesh/process_mesh_data.h"
#include "tools/custom_errors.h"

#include <typeinfo>

namespace DGTD{
  using namespace DG;
  using namespace TD;

template <class Pde, class Basis, class TD_solver>
Dgtd_solver<Pde, Basis, TD_solver>::Dgtd_solver(
    const std::string &_mesh_name,
    const size_t _polynomial_order,
    const double _end_time)
    : mesh_name(_mesh_name),
      polynomial_order(_polynomial_order),
      end_time(_end_time),
      quad_nodes(basis.get_quad_nodes(_polynomial_order)),
      diff_matrix(Elementwise_operations<Basis>(_polynomial_order)
                      .get_diff_matrix()),
      lift_matrix(Elementwise_operations<Basis>(_polynomial_order)
                      .get_lift_matrix()) {

  if (!std::is_same<TD_solver, Low_storage_runge_kutta>::value) {
    throw Not_implemented("Given time-domain solver unknown.");
  }

  //this->time_step = this->get_time_step(end_time);
}


template <class Pde, class Basis, class TD_solver>
arma::mat Dgtd_solver<Pde, Basis, TD_solver>::get_phys_node_coords() const {

  Mesh::Process_mesh_data pmd(mesh_name);
  Geometric_operations go(mesh_name);

  std::vector<size_t> elems;
  for (const auto region: pmd.import_gmsh_regions()) {
    elems = pmd.get_ordered_elems(pmd.get_finite_elems(region));
  }

  arma::mat phys_node_coords(this->quad_nodes.size(), elems.size());
  for (size_t n(0); n<elems.size(); ++n) {
    for (size_t m(0); m<quad_nodes.size(); ++m) {
      phys_node_coords(m,n) = 
        go.convert_ref_to_phys_coord(this->quad_nodes[m], elems[n]);
    }
  }

  return phys_node_coords;
}


template <class Pde, class Basis, class TD_solver>
bool Dgtd_solver<Pde, Basis, TD_solver>::is_field_name_valid(
    const std::string& field_name) const {

  const auto field_names(Pde::get_field_names());
  if (field_names.find(field_name) != field_names.end()) {
    return true;
  } else {
    return false;
  }
}

} // namespace DGTD
