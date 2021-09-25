#include "geometric_operations.h"

namespace DG {

Geometric_operations::Geometric_operations(
    const std::string &mesh_filename)
    : Mesh::Process_mesh_data(mesh_filename) {}
//-------------------------------------------------------------------------
double Geometric_operations::get_geometric_factor(const size_t elem_tag) {

  const double ref_elem_size(2.);
  return ref_elem_size / Mesh::Process_mesh_data::get_elem_size(elem_tag);
}
//-------------------------------------------------------------------------
double Geometric_operations::convert_phys_to_ref_coord(
    const double phys_coord,
    const size_t elem_tag) {

  const auto [left_phys_elem_node, right_phys_elem_node] =
      Mesh::Process_mesh_data::get_elem_coords(elem_tag);
  const double left_ref_elem_node(-1.);
  return left_ref_elem_node + (phys_coord - left_phys_elem_node) *
                                  this->get_geometric_factor(elem_tag);
}
//-------------------------------------------------------------------------
double Geometric_operations::convert_ref_to_phys_coord(
    const double ref_coord,
    const size_t elem_tag) {

  const auto [left_phys_elem_node, right_phys_elem_node] =
      Mesh::Process_mesh_data::get_elem_coords(elem_tag);
  const double left_ref_elem_node(-1.);

  return left_phys_elem_node + (ref_coord - left_ref_elem_node) /
                                   this->get_geometric_factor(elem_tag);
}
//-------------------------------------------------------------------------
} // namespace DG
