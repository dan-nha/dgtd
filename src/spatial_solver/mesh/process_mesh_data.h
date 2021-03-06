#ifndef PROCESS_MESH_DATA_H
#define PROCESS_MESH_DATA_H

#include "import_mesh_data.h"

#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace DG::Mesh {

/**
 * @brief Processing the mesh data imported from a Gmsh file, so that we
 * have access to the data which is relevant to the DG scheme.<br>
 * Note, that in DGTD the strong formulation of a given PDE is solved on
 * size independent unit integrals on each finite elements. Some of the
 * following methods are implemented with regards to the conversion from
 * integrals operating on physical coordinates to those unit (reference)
 * integrals.
 */
class Process_mesh_data : public Import_mesh_data {
public:
  Process_mesh_data(const std::string &mesh_name);

  /// @brief Get finite element tags for a given region
  std::vector<size_t> get_finite_elems(const std::string &region_name);
  std::vector<size_t> get_finite_elems(const size_t region_tag);

  /**
   * @brief Order the regions from left to right. This is useful when
   * operating different solution scheme parameters for different regions
   * and assembling the global solution afterwards.
   */
  std::vector<size_t> get_ordered_regions();

  /**
   * @brief Order a given set of finite elements from left to right. Use
   * this method, to always have the correct order of elements. For
   * example, you can apply this method to a collection of element tags of
   * a certain region. This method is also utilized to figure out the
   * neighboring regions, i.e. the region ordering from left to right.
   */
  std::vector<size_t>
  get_ordered_elems(const std::vector<size_t> &elem_tags);

  /**
   * @brief The element size is needed to convert physical coordinate
   * information to reference coordinate information and vice versa.
   * */
  double get_elem_size(const size_t elem_tag);

  /**
   * @brief Get left and right node coordinate of a 1D finite element.
   * Among others, this method is needed to calculate the element size.
   * */
  std::tuple<double, double> get_elem_coords(const size_t elem_tag);

  /**
   * @brief Get smallest finite element size within the entire mesh. This
   * is needed if I want to calculate the smallest time step for the
   * entire mesh using the CFL critereon.
   */
  double get_min_elem_size();

  /**
   * @brief Get smallest finite element size within a region. This is
   * needed if I want to calculate the smallest time step for a certain
   * region using the CFL critiereon. One could use this method, to adapt
   * the time step to different regional mesh sizes.
   */
  double get_min_elem_size(const std::string &region_name);
  double get_min_elem_size(const size_t region_tag);

private:
  const std::string &mesh_name;
  const size_t dimension;

  double get_coord(const size_t node_tag, const size_t entity_tag);
  std::vector<size_t> get_entity_tags(const size_t region_tag);
  bool is_region(const size_t physical_tag);
  bool is_contour(const size_t physical_tag);
};
} // namespace DG::Mesh

#endif
