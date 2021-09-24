#ifndef PROCESS_MESH_DATA_H
#define PROCESS_MESH_DATA_H

#include "import_mesh_data.h"

#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace DG_solver::Mesh {

/**
 * @brief Processing the mesh data imported from a Gmsh file, so that we
 * have access to the data which is relevant to the DG scheme.
 */
class Process_mesh_data : public Import_mesh_data {
public:
  Process_mesh_data(const std::string &mesh_name);

  /// @brief Get finite element tags for a given region name
  std::vector<size_t> get_finite_elems(const std::string &region_name);
  /// @brief Get finite element tags for a given region
  std::vector<size_t> get_finite_elems(const size_t region_tag);

  /// @brief Order a given set of finite elements from left to right
  std::vector<size_t>
  get_ordered_elems(const std::vector<size_t> elem_tags);

  /// @brief Get left and right node coordinate of a 1D finite element
  std::tuple<double, double> get_elem_coords(const size_t elem_tag);

  /// @brief Get smallest finite element size within the entire mesh
  double get_min_elem_size();

  /// @brief Get smallest finite element size within a region
  double get_min_elem_size(const std::string& region_name);
  /// @brief Get smallest finite element size within a region
  double get_min_elem_size(const size_t region_tag);

  double get_elem_size(const size_t elem_tag);

private:
  const std::string &mesh_name;
  const size_t dimension;

  std::vector<size_t> get_region_tags();
  double get_coord(const size_t node_tag, const size_t entity_tag);
  std::vector<size_t> get_entity_tags(const size_t region_tag);
  bool is_contour(const size_t physical_tag);
  bool is_region(const size_t physical_tag);
};
} // namespace DG_solver::Mesh

#endif
