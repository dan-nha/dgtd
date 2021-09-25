#include "process_mesh_data.h"

namespace DG::Mesh {

Process_mesh_data::Process_mesh_data(const std::string &_mesh_name)
    : Import_mesh_data(_mesh_name), mesh_name(_mesh_name),
      dimension(Import_mesh_data::get_dimension()) {}
//-------------------------------------------------------------------------
std::vector<size_t>
Process_mesh_data::get_finite_elems(const std::string &region_name) {

  auto phys_name_map(Import_mesh_data::import_gmsh_physical_names());
  if (phys_name_map.find(region_name) != phys_name_map.end()) {
    const size_t region_tag(phys_name_map[region_name]);
    if (is_region(region_tag) == false) {
      throw std::invalid_argument(
          std::string() + "__FILE__" + ":" + std::to_string(__LINE__) +
          ": "
          "Invalid region name. Given name belongs to a contour.");
    } else {
      return this->get_finite_elems(region_tag);
    }
  } else {
    throw std::invalid_argument(
        std::string() + "__FILE__" + ":" + std::to_string(__LINE__) +
        ": "
        "Invalid region name.");
  }
}
//----
std::vector<size_t>
Process_mesh_data::get_finite_elems(const size_t region_tag) {

  std::vector<size_t> element_tags;
  for (const auto entity_tag : this->get_entity_tags(region_tag)) {
    for (const auto [element_tag, node_tags] :
         Import_mesh_data::import_gmsh_elements(
             Entity.curve, entity_tag)) {
      element_tags.push_back(element_tag);
    }
  }
  return element_tags;
}
//-------------------------------------------------------------------------
std::vector<size_t>
Process_mesh_data::get_ordered_elems(const std::vector<size_t> elem_tags) {

  std::vector<size_t> ordered_elem_tags;
  std::map<double, size_t> coord_elemtag_map;
  for (const auto elem_tag : elem_tags) {
    auto [left_coord, right_coord] = this->get_elem_coords(elem_tag);
    coord_elemtag_map[left_coord] = elem_tag;
  };

  for (const auto [coord, elem_tag] : coord_elemtag_map) {
    ordered_elem_tags.push_back(elem_tag);
  }
  return ordered_elem_tags;
}
//-------------------------------------------------------------------------
double Process_mesh_data::get_elem_size(const size_t elem_tag) {

  auto [left_coord, right_coord] = this->get_elem_coords(elem_tag);
  return fabs(right_coord - left_coord);
}
//-------------------------------------------------------------------------
std::tuple<double, double>
Process_mesh_data::get_elem_coords(const size_t elem_tag) {

  double left_coord;
  double right_coord;
  for (const auto region_tag : this->get_region_tags()) {
    for (const auto entity_tag : this->get_entity_tags(region_tag)) {
      std::vector<size_t> node_tags(Import_mesh_data::import_gmsh_elements(
          Entity.curve, entity_tag)[elem_tag]);
      if (!node_tags.empty()) {
        try {
          left_coord = this->get_coord(node_tags.front(), entity_tag);
          right_coord = this->get_coord(node_tags.back(), entity_tag);
        } catch (std::invalid_argument &ia) {
          std::cerr << __FILE__ << ":" << __LINE__ << ": "
                    << "Could not find node coordinates for element "
                    << elem_tag << '\n'
                    << ia.what() << std::endl;
        }
      }
    }
  }

  return {left_coord, right_coord};
}
//-------------------------------------------------------------------------
double Process_mesh_data::get_min_elem_size() {

  std::vector<double> elem_sizes;
  for (const auto region_tag : this->get_region_tags()) {
    elem_sizes.push_back(get_min_elem_size(region_tag));
  }

  return *std::min_element(std::begin(elem_sizes), std::end(elem_sizes));
}
//----
double
Process_mesh_data::get_min_elem_size(const std::string &region_name) {

  auto phys_name_map(Import_mesh_data::import_gmsh_physical_names());
  if (phys_name_map.find(region_name) != phys_name_map.end()) {
    const size_t region_tag(phys_name_map[region_name]);
    if (is_region(region_tag) == false) {
      throw std::invalid_argument(
          std::string() + "__FILE__" + ":" + std::to_string(__LINE__) +
          ": "
          "Invalid region name. Given name belongs to a contour.");
    } else {
      return this->get_min_elem_size(region_tag);
    }
  } else {
    throw std::invalid_argument(
        std::string() + "__FILE__" + ":" + std::to_string(__LINE__) +
        ": "
        "Invalid region name.");
  }
}
//----
double Process_mesh_data::get_min_elem_size(const size_t region_tag) {

  std::vector elem_tags(this->get_finite_elems(region_tag));

  std::vector<double> elem_sizes;
  for (const auto elem_tag : elem_tags) {
    elem_sizes.push_back(this->get_elem_size(elem_tag));
  }

  return *std::min_element(std::begin(elem_sizes), std::end(elem_sizes));
}
//-------------------------------------------------------------------------
std::vector<size_t> Process_mesh_data::get_region_tags() {

  const size_t region_entity(Mesh::physgroup_entity_table(
      this->dimension)[Physical_group.region]);

  std::vector<size_t> physical_tags;
  for (const auto [physical_tag, entity_dim] :
       Import_mesh_data::import_gmsh_physical_groups()) {
    if (entity_dim == region_entity) {
      physical_tags.push_back(physical_tag);
    }
  }

  return physical_tags;
}
//-------------------------------------------------------------------------
double Process_mesh_data::get_coord(
    const size_t node_tag,
    const size_t entity_tag) {
  auto curve_node(Import_mesh_data::import_gmsh_nodes());
  if (curve_node.find(node_tag) != curve_node.end()) {
    return curve_node[node_tag].front();
  } else {
    throw std::invalid_argument(
        std::string() + __FILE__ + ":" + std::to_string(__LINE__) +
        ": "
        "Did not find coordinate for given node and entity.");
  }
}
//-------------------------------------------------------------------------
std::vector<size_t>
Process_mesh_data::get_entity_tags(const size_t region_tag) {

  std::vector<size_t> entity_tags;
  for (const auto [entity_tag, physical_tags] :
       Import_mesh_data::import_gmsh_entities(Entity.curve)) {
    for (const auto physical_tag : physical_tags) {
      if (physical_tag == region_tag) {
        entity_tags.push_back(entity_tag);
      }
    }
  }
  return entity_tags;
}
//-------------------------------------------------------------------------
bool Process_mesh_data::is_contour(const size_t physical_tag) {

  const size_t entity(
      Import_mesh_data::import_gmsh_physical_groups()[physical_tag]);

  if (Mesh::entity_physgroup_table(this->dimension)[entity] ==
      Physical_group.contour) {
    return true;
  } else {
    return false;
  }
}
//-------------------------------------------------------------------------
bool Process_mesh_data::is_region(const size_t physical_tag) {

  const size_t entity(
      Import_mesh_data::import_gmsh_physical_groups()[physical_tag]);

  if (Mesh::entity_physgroup_table(this->dimension)[entity] ==
      Physical_group.region) {
    return true;
  } else {
    return false;
  }
}

} // namespace DG::Mesh
