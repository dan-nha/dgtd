#ifndef IMPORT_MESH_DATA_H
#define IMPORT_MESH_DATA_H

#include "mesh.h"

#include <armadillo>
#include <fstream>
#include <tuple>
#include <vector>

namespace DG_solver::Mesh {

class Import_mesh_data {

public:
  Import_mesh_data(const std::string &filename);

  /**
   * @brief Extract the overall dimension of the mesh from the section
   * '$PhysicalNames' in the mesh file. The highest dimension found under
   * '$PhysicalNames' is equal to the overall dimension of the mesh.
   *
   * @return Dimension of the mesh
   */
  size_t get_dimension() const;

  /// @brief Import the number of mesh nodes (not quadrature nodes!)
  size_t import_number_of_nodes() const;

protected:
  /** 
   * @brief Import Gmsh "$PhysicalNames" map: "physicalTag" (key) to
   * "dimension" (value)
   */
  std::map<size_t, size_t> import_gmsh_physical_groups();

  /**
   * @brief For a given entity (point, curve, surface, volume) a map of the
   * entity tags within each entity type and the physical tags associated
   * with each entity tag is generated
   *
   * @param[in] entity_type The type of entity, i.e. Entity.point,
   * Entity.curve, Entity.surface, or Entity.volume according to gmsh
   *
   * @return Map of gmsh "entityTag" (key) to a vector of "physicalTags"
   * (value)
   */
  std::map<size_t, std::vector<size_t>>
  import_gmsh_entities(const size_t entity_type);

  /**
   * @brief Import the mesh nodes for a given entity, not to be confused
   * with the quadrature nodes in DGTD. The entity is characterized by its
   * type and tag number.
   *
   * @param[in] entity_type The type of entity, i.e. Entity.point,
   * Entity.curve, Entity.surface, or Entity.volume according to gmsh
   * @param[in] entity_tag The gmsh entity tag, i.e. the number
   * characterizing the entity
   *
   * @return Map of gmsh "nodeTag" (key) to the node's coordinates (value).
   */
  std::map<size_t, arma::vec>
  import_gmsh_nodes(const size_t entity_type, const size_t entity_tag);

  /**
   * @brief Import gmsh elements for a certain entity type and tag. Note,
   * that gmsh elements are not elements in the sense of elements in finite
   * elements but rather a collection of points and curves in 1D, curves
   * and sufaces in 2D, or surfaces and volumes in 3D which form a certain
   * entity.
   *
   * @param[in] entity_type The type of entity, i.e. Entity.point,
   * Entity.curve, Entity.surface, or Entity.volume
   * @param[in] entity_tag The entity tag, i.e. the number characterizing
   * the entity
   *
   * @return Map of gmsh "elementTag" (key) to the associated "nodeTags"
   * (value), i.e. the tags of those nodes which form the element
   */
  std::map<size_t, std::vector<size_t>>
  import_gmsh_elements(const size_t entity_type, const size_t entity_tag);

  /**
   * @brief Import gmsh element type for a certain entity type and tag.
   * Note, that gmsh elements are not elements in the sense of elements in
   * finite elements but rather a collection of points and curves in 1D,
   * curves and sufaces in 2D, or surfaces and volumes in 3D which form a
   * certain entity.
   *
   * @param[in] entity_type The type of entity, i.e. Entity.point,
   * Entity.curve, Entity.surface, or Entity.volume
   * @param[in] entity_tag The entity tag, i.e. the number characterizing
   * the entity
   *
   * @return Element type as integer
   */
  size_t import_gmsh_element_type(
      const size_t entity_type,
      const size_t entity_tag);

private:
  /**
   * @brief Mesh section information is given by: the mesh specifier name
   * from the specifier_list in mesh.h the position in the mesh, where the
   * specifier is found the line number in the mesh, where the specifier is
   * found
   */
  typedef std::tuple<std::string, size_t, size_t> mesh_section_info;

  /**
   * @brief Extract the mesh section information from the mesh. The
   * sections are given according to the specifier_list in mesh.h. More
   * specifically, the sections are specified by "$MeshFormat",
   * "$PhysicalNames", "$Entities", "$Nodes", and "$Elements"
   *
   * @return Mesh specifier name, position of this specifier in the mesh,
   * line number for each mesh section
   */
  std::vector<mesh_section_info> get_mesh_section_info() const;

  /**
   * @brief Go to a specific section in the mesh, given by one of the
   * specfiers "$MeshFormat", "$PhysicalNames", "$Entities", "$Nodes", or
   * "$Elements"
   *
   * @param[in] specifier Specifier of the mesh section
   * @param[in] mesh_file Mesh file to read from
   */
  void goto_mesh_section(
      const std::string &specifier,
      std::ifstream &mesh_file) const;

  /**
   * @brief Go to a specific section in the mesh, given by one of the
   * specfiers "$MeshFormat", "$PhysicalNames", "$Entities", "$Nodes", or
   * "$Elements", while incrementing a given line number
   *
   * @param[in] specifier Specifier of the mesh section
   * @param[in] mesh_file Mesh file to read from
   * @param[in] line_number Line number to increment by going to the mesh
   * section
   */
  void goto_mesh_section(
      const std::string &specifier,
      std::ifstream &mesh_file,
      size_t &line_number) const;

  const std::string mesh_name;
  const std::string mesh_error;

protected:
  const size_t dimension;
};
} // namespace DG_solver::Mesh

#endif
