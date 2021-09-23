#include "check_mesh.h"
#include "../../tools/import.h"
#include "../../tools/custom_errors.h"

#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/log/trivial.hpp>

namespace Mesh {
Check_mesh::Check_mesh(const std::string &filename)
    : mesh_name(filename),
    begin_specifier({
      {"$MeshFormat", false},
      {"$PhysicalNames", false},
      {"$Entities", false},
      {"$Nodes", false},
      {"$Elements", false}}),
    end_specifier({
      {"$EndMeshFormat", false},
      {"$EndPhysicalNames", false},
      {"$EndEntities", false},
      {"$EndNodes", false},
      {"$EndElements", false}})
{
  this->check_mesh_existence();
  this->check_empty_lines();
  this->check_specifiers();
  this->check_mesh_format();
}
//-------------------------------------------------------------------------
void Check_mesh::check_mesh_existence() {

  std::ifstream mesh_file(this->mesh_name.c_str());
  boost::filesystem::path cur_file_path(mesh_name);
  std::string file_extension(boost::filesystem::extension(cur_file_path));

  if (file_extension != ".msh") {
    throw Mesh_error(
        "Mesh file seems not to be of type '.msh'", this->mesh_name);
  }

  if (!mesh_file.good()) {
    throw std::ifstream::failure(
        "Error loading file '" + this->mesh_name + "'.");
  }

  mesh_file.close();
}
//-------------------------------------------------------------------------
void Check_mesh::check_empty_lines() {
  std::ifstream mesh_file(this->mesh_name.c_str());

  if (mesh_file.peek() == std::ifstream::traits_type::eof()) {
    throw Mesh_error("File seems to be empty.", this->mesh_name);
  }

  std::string line;
  for (size_t line_number(1); std::getline(mesh_file, line);
       ++line_number) {
    if (line.empty()) {
      throw Mesh_error(
          "Detected empty line in mesh file. (line " +
              std::to_string(line_number) + ")",
          this->mesh_name);
    }
  }
}
//-------------------------------------------------------------------------
void Check_mesh::check_specifiers() {

  std::ifstream mesh_file(this->mesh_name.c_str());
 
  // Look for begin and end specifiers in mesh file
  // ----------------------------------------------
  size_t line_number(1);
  for (std::string line; std::getline(mesh_file, line); ++line_number) {
    const std::string first_word(Import::get_entry<std::string>(line));

    auto begin_iter = this->begin_specifier.begin();
    auto end_iter = this->end_specifier.begin();
    for (; (begin_iter != this->begin_specifier.end()) &&
           (end_iter != this->end_specifier.end());
         ++begin_iter, ++end_iter) {

      // Look for begin specifier
      if (first_word == (*begin_iter).first) {

        (*begin_iter).second = true;
        size_t pos = mesh_file.tellg();

        // Check whether the next line is already the end specifer and, if
        // so,throw error
        std::string next_line;
        std::getline(mesh_file, next_line);
        if (next_line == (*end_iter).first) {
          throw Mesh_error(
              "Missing content between " + (*begin_iter).first +
                  " specifiers.",
              this->mesh_name);
        }
        mesh_file.seekg(pos);
      }

      // Look for end specifier
      if (first_word == (*end_iter).first)
        (*end_iter).second = true;
    }
  }

  this->check_unfound_specifiers();
}
//-------------------------------------------------------------------------
void Check_mesh::check_mesh_format() {

  /// Checks the mesh version tag
  /// ==============================
  std::ifstream mesh_file(this->mesh_name.c_str());

  for (std::string line; std::getline(mesh_file, line);) {
    if (line == "$MeshFormat")
      break;
  }

  /// Parse and check mesh version
  const double mesh_version(
      Import::get_next_line_entry<double>(mesh_file));

  if (mesh_version != 4.1) {
    throw Mesh_error(
        "Mesh format '" + std::to_string(mesh_version) +
            "' is not supported. "
            "Only meshes of format 4.1 produced by gmsh are currently "
            "supported. "
            "Try to export your gmsh mesh by using 'gmsh -format msh41'.",
        this->mesh_name);
  }
}
//-------------------------------------------------------------------------
void Check_mesh::check_unfound_specifiers() {
  const std::string unfound_specifiers(this->get_unfound_specifiers());
  if (!(unfound_specifiers.empty())) {
    throw Mesh_error(
        "The following specifiers are missing in the mesh file:\n" +
            unfound_specifiers,
        this->mesh_name);
  }
}
//-------------------------------------------------------------------------
std::string Check_mesh::get_unfound_specifiers() const {
  std::string unfound_specifiers;

  for (auto &unfound : this->begin_specifier) {
    if (!unfound.second)
      unfound_specifiers.append(unfound.first + '\n');
  }

  for (auto &unfound : this->end_specifier) {
    if (!unfound.second)
      unfound_specifiers.append(unfound.first + '\n');
  }
  
  return unfound_specifiers;
}
} // namespace Mesh
