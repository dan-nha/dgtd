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
      const double end_time);

  arma::mat get_phys_node_coords() const;

  bool is_field_name_valid(const std::string &field_name) const;

private:
  Basis basis;
  const std::string mesh_name;
  const size_t polynomial_order;
  const double end_time;
  std::vector<double> quad_nodes;
  const arma::mat diff_matrix;
  const arma::mat lift_matrix;
};
} // namespace DGTD

#include "dgtd_solver.tpp"

#endif
