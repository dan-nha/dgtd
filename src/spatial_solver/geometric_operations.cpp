#include "geometric_operations.h"

#include <cmath>
#include <algorithm>

namespace DG {

Geometric_operations::Geometric_operations(
    Mesh::Process_mesh_data &_processed_mesh)
    : processed_mesh{_processed_mesh} {}
//-------------------------------------------------------------------------
double Geometric_operations::get_geometric_factor(const size_t elem_tag) {

  const double ref_elem_size{2.};
  return ref_elem_size / processed_mesh.get_elem_size(elem_tag);
}
//-------------------------------------------------------------------------
double Geometric_operations::convert_phys_to_ref_coord(
    const double phys_coord,
    const size_t elem_tag) {

  const auto [left_phys_elem_node, right_phys_elem_node] =
    processed_mesh.get_elem_coords(elem_tag);
  const double left_ref_elem_node{-1.};
  return left_ref_elem_node + (phys_coord - left_phys_elem_node) *
                                  this->get_geometric_factor(elem_tag);
}
//-------------------------------------------------------------------------
double Geometric_operations::convert_ref_to_phys_coord(
    const double ref_coord,
    const size_t elem_tag) {

  const auto [left_phys_elem_node, right_phys_elem_node] =
    processed_mesh.get_elem_coords(elem_tag);
  const double left_ref_elem_node{-1.};

  return left_phys_elem_node + (ref_coord - left_ref_elem_node) /
                                   this->get_geometric_factor(elem_tag);
}
//-------------------------------------------------------------------------
double Geometric_operations::get_min_node_dist(
    const std::vector<double> quad_nodes) {

  const double ref_elem_size{2.};
  const double min_elem_size(processed_mesh.get_min_elem_size());
  return this->get_min_distance(quad_nodes) * min_elem_size / ref_elem_size;
}
//-------------------------------------------------------------------------
double Geometric_operations::get_min_node_dist(
    const std::vector<double> quad_nodes,
    const size_t region_tag) {

  const double ref_elem_size{2.};
  const double min_elem_size(processed_mesh.get_min_elem_size(region_tag));
  return this->get_min_distance(quad_nodes) * min_elem_size / ref_elem_size;
}
//-------------------------------------------------------------------------
double Geometric_operations::get_min_distance(
    const std::vector<double> coords) const {
  
  std::vector<double> distances;
  for (size_t i{0}; i<coords.size()-1; ++i) {
    distances.push_back(fabs(coords[i+1]-coords[i]));
  }
  return *std::min_element(distances.begin(), distances.end());
}

} // namespace DG
