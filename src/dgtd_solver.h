#ifndef DGTD_SOLVER_H
#define DGTD_SOLVER_H

#include "tools/input.h"
#include "spatial_solver/mesh/process_mesh_data.h"
#include "spatial_solver/geometric_operations.h"

#include <armadillo>
#include <string>

namespace DGTD {

template <class Pde, class Basis, class TD_solver> class Dgtd_solver {
public:
  Dgtd_solver(DG::Mesh::Process_mesh_data &processed_mesh, 
      const Input &input);

  arma::mat get_solution(Pde pde);

  arma::mat get_phys_node_coords();

  // Get geometric factors for all elements
  std::vector<double> get_geometric_factors();

  double get_time_step();

  double get_min_node_dist();

  bool is_field_name_valid(const std::string &field_name) const;

private:
  DG::Mesh::Process_mesh_data &processed_mesh;
  const Input &input;
  DG::Geometric_operations geo;
  const Basis basis;
  const size_t polynomial_order;
  const std::vector<double> quad_nodes;
  const double end_time;
  const double dt_factor;
  const double time_step;
  const arma::mat diff_matrix;
  const arma::mat lift_matrix;
};
} // namespace DGTD

#include "dgtd_solver.tpp"

#endif
