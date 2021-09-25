#ifndef GEOMETRY_OPERATIONS_H
#define GEOMETRY_OPERATIONS_H

#include "mesh/process_mesh_data.h"

#include <string>
#include <vector>

namespace DG_solver {
/**
 * @brief In DGTD the strong formulation of a given PDE is solved on
 * size independent unit integrals on each finite elements. Some of the
 * following methods are implemented with regards to the conversion from
 * integrals operating on physical coordinates to those unit (reference)
 * integrals.
 */
class Geometric_operations : public Mesh::Process_mesh_data {
public:
  Geometric_operations(const std::string &mesh_filename);

  /**
   * @brief Conversion factor (Jacobian) between reference unit integral
   * and integral on physical node coordinates of a finite element
   */
  double get_geometric_factor(const size_t elem_tag);

  /**
   * @brief Coordinate conversion from physical coordinates of an element
   * to references coordinates of a unit element with a length of 2
   */
  double convert_phys_to_ref_coord(
      const double phys_coord,
      const size_t elem_tag);

  /// @brief Back-transform of convert_phys_to_ref_coord(...)
  double convert_ref_to_phys_coord(
      const double ref_coord, 
      const size_t elem_tag);
};
} // namespace DG_solver

#endif
