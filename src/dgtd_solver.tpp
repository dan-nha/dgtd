#include "dgtd_solver.h"
#include "spatial_solver/Elementwise_operations.h"
#include "spatial_solver/Geometric_operations.h"
#inlcude "temporal_solver/low_storage_runge_kutta.h"

namespace DG::TD {
template <class Pde, class Basis, class TD_solver>
Dgtd_solver<Pde, Basis, TD_solver>::Dgtd_solver(
    const std::string &mesh_name,
    const size_t _polynomial_order,
    const double _end_time)
    : polynomial_order(_polynomial_order),
      diff_matrix(Elementwise_operations<Basis>(polynomial_order_)
                      .get_differentiation_matrix()),
      lift_matrix(Elementwise_operations<Basis>(polynomial_order_)
                      .get_lift_matrix()),
      end_time(_end_time) {

  if (!std::is_same<TD_solver, Temporal_solver::Low_storage_runge_kutta>::
          value) {
    throw Not_implemented("Given time-domain solver unknown.");
  }

  this->time_step = this->get_time_step(end_time);
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

} // namespace DG::TD
