#ifndef DGTD_SOLVER_H
#define DGTD_SOLVER_H

#include <string>
#include <armadillo>

namespace DGTD{
template <class Pde, class Basis, class TD_solver> class Dgtd_solver {
public:
  Dgtd_solver(
      const std::string &mesh_name,
      const size_t polynomial_order,
      const double end_time,
      const double dt_factor,
      const double upwind_param);

  arma::mat get_solution(
      Pde pde,
      const size_t runge_kutta_order, 
      const size_t runge_kutta_stages) const;
  arma::mat get_phys_node_coords() const;
  std::vector<double> get_geometric_factors() const;
  double get_time_step();
  double get_min_node_dist();
  bool is_field_name_valid(const std::string &field_name) const;
  void store_results() const;

private:
  Basis basis;
  const std::string mesh_name;
  const size_t polynomial_order;
  std::vector<double> quad_nodes;
  const double end_time;
  const double dt_factor;
  const double time_step;
  const double upwind_param;
  const arma::mat diff_matrix;
  const arma::mat lift_matrix;
};
} // namespace DGTD

#include "dgtd_solver.tpp"

#endif
