#include "mesh_geometry.h"
#include "import_mesh_data.h"
#include "../../tools/import.h"
#include "../../tools/custom_errors.h"

#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <climits>
#include <sstream>

namespace DG_solver::Mesh {

Import_mesh_data::Import_mesh_data(const std::string &filename)
    : mesh_name(filename), dimension(this->get_dimension()) {}
//-----------------------------------------------------------------------
size_t Import_mesh_data::get_dimension() const {

  std::ifstream mesh_file(mesh_name.c_str());
  this->goto_mesh_section("$PhysicalNames", mesh_file);

  Import::skip_lines(mesh_file, 1);

  std::vector<size_t> phys_group_dims;
  for (std::string line; std::getline(mesh_file, line);) {

    std::string line_entry(Import::get_entry<std::string>(line));

    if (line_entry == "$EndPhysicalNames")
      break;

    phys_group_dims.push_back(boost::lexical_cast<size_t>(line_entry));

    if (phys_group_dims.back() > 3) {
      throw Mesh_error(
          std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
              ": "
              "Invalid dimension.",
          this->mesh_name);
    } 
  }

  return *std::max_element(phys_group_dims.begin(), phys_group_dims.end());
}
//-----------------------------------------------------------------------
size_t Import_mesh_data::import_number_of_nodes() const {

  std::ifstream mesh_file(mesh_name.c_str());
  Import_mesh_data::goto_mesh_section("$Nodes", mesh_file);

  return Import::get_next_line_entry<size_t>(mesh_file, 2);
}
//-----------------------------------------------------------------------
std::map<size_t, size_t> Import_mesh_data::import_gmsh_physical_groups() {

  std::ifstream mesh_file(this->mesh_name.c_str());
  this->goto_mesh_section("$PhysicalNames", mesh_file);

  const size_t phys_group_number(
      Import::get_next_line_entry<size_t>(mesh_file));

  std::map<size_t, size_t> physical_group_info;

  for (size_t phys_group_counter = 0;
       phys_group_counter < phys_group_number;
       ++phys_group_counter) {

    std::string line;
    std::getline(mesh_file, line);
    size_t phys_group_dim(Import::get_entry<size_t>(line));
    size_t phys_group_tag(Import::get_entry<size_t>(line, 2));

    physical_group_info[phys_group_tag] = phys_group_dim;
  }

  if (Import::get_next_line_entry<std::string>(mesh_file) !=
      "$EndPhysicalNames") {
    BOOST_LOG_TRIVIAL(warning)
        << __FILE__ << ":" << __LINE__
        << "Ignoring excess entries in '$PhysicalNames'." << std::endl;
  }

  return physical_group_info;
}
//-----------------------------------------------------------------------
std::map<size_t, std::vector<size_t>>
Import_mesh_data::import_gmsh_entities(const size_t entity_type) {

  std::ifstream mesh_file(this->mesh_name.c_str());
  size_t line_number; // Only used for warning message below
  this->goto_mesh_section("$Entities", mesh_file, line_number);

  // Import number of point, curve, surface, and volume entities in the
  // mesh
  std::vector<size_t> number_of_each_entity(
      Import::get_next_line_entries<size_t>(mesh_file, line_number));
  size_t &num_entities(number_of_each_entity[entity_type]);

  // Set the line and position in the line from which to read the physical
  // tags from depending on the entity type
  size_t phys_tag_pos;
  if (entity_type == Entity.point) {
    phys_tag_pos = 5;
  } else {
    phys_tag_pos = 8;

    size_t number_of_previous_entities(0);
    for (size_t prev_entity_type(0); prev_entity_type < entity_type;
         ++prev_entity_type) {
      number_of_previous_entities +=
          number_of_each_entity[prev_entity_type];
    }
    Import::skip_lines(mesh_file, number_of_previous_entities);
    line_number += number_of_previous_entities;
  }

  // Import entity tag and associated physical tags for this entity
  std::map<size_t, std::vector<size_t>> entity_info;
  std::string line;
  for (size_t entity_tag_idx(0);
       entity_tag_idx < num_entities && getline(mesh_file, line);
       ++entity_tag_idx, ++line_number) {

    size_t entity_tag(Import::get_entry<size_t>(line, 1));

    // BUG in gmsh: in 2D and 3D the 8-th line entry is not a size_t but a
    // double contrary to the documentation
    // (https://gmsh.info/doc/texinfo/gmsh.html). Hence, a double is
    // imported and casted onto a size_t
    size_t num_phys_tags(Import::get_entry<double>(line, phys_tag_pos));

    if (num_phys_tags > 1 && entity_type == this->dimension) {
      throw Mesh_error(
          std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
              ": Region entity with more than one physical tag detected. "
              "(line " +
              std::to_string(line_number) + " in mesh)",
          this->mesh_name);
    };

    std::vector<size_t> physical_tags;
    if (num_phys_tags == 0 && entity_type > Entity.point) {
      physical_tags.push_back(
          // BUG in gmsh: same bug as above
          Import::get_entry<double>(line, phys_tag_pos + 1));
    } else {
      for (size_t phys_tag_idx(1); phys_tag_idx <= num_phys_tags;
           ++phys_tag_idx) {
        physical_tags.push_back(
            // BUG in gmsh: same bug as above
            Import::get_entry<double>(line, phys_tag_pos + phys_tag_idx));
      }
    }

    entity_info[entity_tag] = physical_tags;
  }

  return entity_info;
}
//-----------------------------------------------------------------------
std::map<size_t, arma::vec> Import_mesh_data::import_gmsh_nodes(
    const size_t entity_type,
    const size_t entity_tag) {

  std::ifstream mesh_file(this->mesh_name.c_str());
  this->goto_mesh_section("$Nodes", mesh_file);

  size_t num_entity_blocks(Import::get_next_line_entry<size_t>(mesh_file));

  std::map<size_t, arma::vec> node_map;
  for (size_t entity_block(1); entity_block <= num_entity_blocks;
       ++entity_block) {

    std::vector<size_t> block_info(
        Import::get_next_line_entries<size_t>(mesh_file));
    size_t block_entity_dim, block_entity_tag, parametric_flag,
        block_num_nodes;
    if (block_info.size() == 4 && block_info[0] <= this->dimension &&
        block_info[2] <= 1) {
      block_entity_dim = block_info[0];
      block_entity_tag = block_info[1];
      parametric_flag = block_info[2];
      block_num_nodes = block_info[3];

      if (parametric_flag == 1) {
        throw Mesh_error(
            std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
                ": Non-trivial parametric flag in $Nodes section entity "
                "block #" +
                std::to_string(entity_block) +
                " detected. Mesh parametrization is currently not "
                "supported in miniDGTD. Regenerate mesh.",
            this->mesh_name);
      }
    }

    if (block_entity_dim == entity_type &&
        block_entity_tag == entity_tag) {

      if (block_num_nodes != 0) {
        // Get the gmsh "nodeTags"
        std::vector<size_t> node_tags;
        for (size_t node(0); node < block_num_nodes; ++node) {
          // BUG in gmsh:
          // Reading in double instead of size_t and then perform a static
          // cast to size_t due to bug in 2D and 3D mesh file produced by
          // gmsh
          node_tags.push_back(static_cast<size_t>(
              Import::get_next_line_entry<double>(mesh_file)));
        }

        // Get the node coordinates
        for (const size_t &node_tag : node_tags) {
          arma::vec coords(
              Import::get_next_line_entries<double>(mesh_file));
          node_map[node_tag] = coords;
        }

      } else {
        std::string error_message(
            std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
            ": No nodes for entity type " + std::to_string(entity_type) +
            " and entity tag " + std::to_string(entity_tag) + ".");
        throw Mesh_error(error_message, this->mesh_name);
      }
    } else {
      Import::skip_lines(mesh_file, 2 * block_num_nodes);
    }
  }

  if (node_map.empty()) {
    throw Mesh_error(
        std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
        ": Unknown entity type " + std::to_string(entity_type) +
        " or tag " + std::to_string(entity_tag) + ".",
        this->mesh_name);
  }

  return node_map;
}
//-----------------------------------------------------------------------
std::map<size_t, std::vector<size_t>>
Import_mesh_data::import_gmsh_elements(
    const size_t entity_type,
    const size_t entity_tag) {

  std::ifstream mesh_file(this->mesh_name.c_str());
  this->goto_mesh_section("$Elements", mesh_file);

  std::map<size_t, std::vector<size_t>> element_info;

  size_t num_entity_blocks(Import::get_next_line_entry<size_t>(mesh_file));

  for (size_t entity_block(0); entity_block < num_entity_blocks;
       ++entity_block) {
    std::vector<size_t> line_entries(
        Import::get_next_line_entries<size_t>(mesh_file));

    if (line_entries.size() == 4 && line_entries[0] <= Entity.volume &&
        line_entries[2] <=
            Element.point // Might need adjustment, if new element types
                          // are introduced in mesh.h
    ) {

      if (line_entries[0] == entity_type &&
          line_entries[1] == entity_tag) {

        size_t num_elements_in_block(line_entries[3]);
        for (size_t element(0); element < num_elements_in_block;
             ++element) {
          std::vector<size_t> element_and_nodes(
              Import::get_next_line_entries<size_t>(mesh_file));
          size_t &element_tag(element_and_nodes[0]);
          std::vector<size_t> node_tags(
              element_and_nodes.begin() + 1, element_and_nodes.end());

          element_info[element_tag] = node_tags;
        }
        break;
      } else {
        Import::skip_lines(mesh_file, line_entries[3]);
      }
    }
  }

  return element_info;
}
//-----------------------------------------------------------------------
size_t Import_mesh_data::import_gmsh_element_type(
    const size_t entity_type,
    const size_t entity_tag) {
  std::ifstream mesh_file(this->mesh_name.c_str());
  this->goto_mesh_section("$Elements", mesh_file);

  size_t element_type;

  size_t num_entity_blocks(Import::get_next_line_entry<double>(mesh_file));

  for (size_t entity_block(0); entity_block < num_entity_blocks;
       ++entity_block) {
    std::vector<size_t> line_entries(
        Import::get_next_line_entries<size_t>(mesh_file));

    if (line_entries.size() == 4 && line_entries[0] == entity_type &&
        line_entries[1] == entity_tag &&
        (line_entries[2] == Element.line_2nodes ||
         line_entries[2] == Element.line_3nodes_2nd_order ||
         line_entries[2] == Element.triangle_3nodes ||
         line_entries[2] == Element.triangle_6nodes_2nd_order ||
         line_entries[2] == Element.tetrahedron_4nodes ||
         line_entries[2] == Element.tetrahedron_10nodes_2nd_order)) {
      element_type = line_entries[2];
      break;
    } else {
      Import::skip_lines(mesh_file, line_entries[3]);
    }
  }
  if (element_type == INT_MAX) {
    throw Mesh_error(
        std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
            "Gmsh element not found.",
        this->mesh_name);
  }

  return element_type;
}
//-----------------------------------------------------------------------
std::vector<Import_mesh_data::mesh_section_info>
Import_mesh_data::get_mesh_section_info() const {

  std::vector<mesh_section_info> mesh_section_info;

  std::fstream mesh_file(this->mesh_name.c_str());
  std::string line;

  for (size_t line_number(1); std::getline(mesh_file, line);
       ++line_number) {

    const std::string first_word(Import::get_entry<std::string>(line));

    for (const auto &specifier : Mesh::specifier_list) {
      if (first_word == specifier) {
        mesh_section_info.push_back(
            {specifier, mesh_file.tellg(), line_number});
      }
    }
  }

  return mesh_section_info;
}
//-----------------------------------------------------------------------
void Import_mesh_data::goto_mesh_section(
    const std::string &specifier,
    std::ifstream &mesh_file) const {

  std::vector<mesh_section_info> section_info(
      this->get_mesh_section_info());
  std::string cur_specifier;
  size_t cur_mesh_pos;

  for (const auto &list_entry : section_info) {
    std::tie(cur_specifier, cur_mesh_pos, std::ignore) = list_entry;
    if (cur_specifier == specifier) {
      mesh_file.seekg(cur_mesh_pos);
      break;
    }
  }
}
//-----------------------------------------
void Import_mesh_data::goto_mesh_section(
    const std::string &specifier,
    std::ifstream &mesh_file,
    size_t &line_number) const {

  std::vector<mesh_section_info> section_info(
      this->get_mesh_section_info());
  std::string cur_specifier;
  size_t cur_mesh_pos;
  size_t cur_line_num;

  for (const auto &list_entry : section_info) {
    std::tie(cur_specifier, cur_mesh_pos, cur_line_num) = list_entry;
    if (cur_specifier == specifier) {
      mesh_file.seekg(cur_mesh_pos);
      line_number = cur_line_num;
      break;
    }
  }
}

} // namespace DG_solver::Mesh
