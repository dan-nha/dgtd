#ifndef DGTD_SOLVER_H
#define DGTD_SOLVER_H
#define BOOST_BIND_GLOBAL_PLACEHOLDERS

#include "tools/input.h"
#include "spatial_solver/mesh/process_mesh_data.h"
#include "spatial_solver/geometric_operations.h"
#include "spatial_solver/elementwise_operations.h"

#include <armadillo>
#include <string>
#include <map>

namespace DGTD {
using namespace DG;

template <class Pde, class Basis, class TD_solver> class Dgtd_solver {
public:
  Dgtd_solver(Mesh::Process_mesh_data &processed_mesh, 
      const Input &input);

  arma::mat get_solution(Pde &pde);

  void initialize_dg_scheme(
      Pde &pde,
      std::map<size_t, arma::mat> &region_field,
      std::map<size_t, arma::mat> &region_phys_node_coords,
      std::map<size_t, std::vector<double>> &region_geo_factors);

  arma::mat evolve_dg_scheme(
      Pde &pde,
      TD_solver &lsrk,
      const arma::mat &fields,
      const double time,
      const std::vector<double> &geo_factors) const;

  arma::mat assemble_global_solution(
      std::map<size_t, arma::mat> &region_field,
      const size_t num_nodes,
      const size_t num_elems);

  arma::mat get_phys_node_coords();
  arma::mat get_phys_node_coords(const size_t region);

  std::vector<double> get_geometric_factors();
  std::vector<double> get_geometric_factors(const size_t region);

  double get_time_step();

  double get_min_node_dist();

  bool is_field_name_valid(const std::string &field_name) const;

private:
  Mesh::Process_mesh_data &processed_mesh;
  const Input &input;
  Geometric_operations geop;
  Elementwise_operations<Basis> eop;
  const Basis basis;
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
