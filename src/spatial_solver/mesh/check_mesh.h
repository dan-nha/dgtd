#ifndef CHECK_MESH_H
#define CHECK_MESH_H

#include "mesh.h"

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

namespace Mesh {

/**
 * @brief Class to check the mesh file in general, more detailed checks can
 * be found in the classes Import_mesh_data and Process_mesh_data
 */
class Check_mesh {

public:
  Check_mesh(const std::string &filename);

  void check_mesh_existence();

  void check_empty_lines();

  void check_mesh_format();

  /**
   * @brief In Gmsh the mesh file is ordered in sections, which are
   * enclosed by so-called specifiers. Hence, the check for these
   * specifiers.
   */
  void check_specifiers();

private:
  /**
   * @brief Use vector of pairs to build a map of specific order. I did not
   * std::map for that matter as the internal ordering differs from the
   * aspired one.
   */
  typedef std::vector<std::pair<std::string, bool>> specifier_table;
  specifier_table begin_specifier;
  specifier_table end_specifier;

  const std::string mesh_name;

  void check_unfound_specifiers();
  std::string get_unfound_specifiers() const;
};

} // namespace Mesh
#endif
