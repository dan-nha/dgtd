#include "../../../src/spatial_solver/mesh/check_mesh.h"
#include "../../../src/tools/custom_errors.h"

#include <boost/test/unit_test.hpp>
#include <string>

namespace Mesh {

/**
 * In the following I check a given mesh for the existence of the most
 * crucial parts of the mesh, which are needed for the DGTD code to run. 
 * I will not check the content of the mesh itself, as I assume that the
 * output of Gmsh is tested. Though minor bugs may occur, they most
 * probably will not affect the DGTD code.
 */
BOOST_AUTO_TEST_SUITE(check_mesh);

const std::string root_dir(DGTD_ROOT);
const std::string
    mesh_dir("/test/spatial_solver/mesh/test_meshes/fail_meshes/");

std::string error_msg;
bool expected_msg(const Mesh_error &ex) {
  BOOST_CHECK_EQUAL(ex.what(), error_msg);
  return true;
};

BOOST_AUTO_TEST_CASE(check_extension) {
  error_msg = "[test.mesh] Mesh file seems not to be of type '.msh'";
  BOOST_CHECK_EXCEPTION(Check_mesh("test.mesh"), Mesh_error, expected_msg);
}

BOOST_AUTO_TEST_CASE(content_existence) {
  const std::string mesh(root_dir + mesh_dir + "empty.msh");
  error_msg = "[" + mesh + "] File seems to be empty.";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(std::string(mesh)), Mesh_error, expected_msg);
}

BOOST_AUTO_TEST_CASE(empty_lines) {
  const std::string mesh(root_dir + mesh_dir + "empty_lines.msh");
  error_msg = "[" + mesh + "] Detected empty line in mesh file. (line 1)";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(std::string(mesh)), Mesh_error, expected_msg);
}

BOOST_AUTO_TEST_CASE(meshformat_specifiers) {
  std::string mesh(
      root_dir + mesh_dir + "missing_meshformat_begin_spec.msh");
  error_msg = "[" + mesh +
              "] The following specifiers are "
              "missing in the mesh file:\n$MeshFormat\n";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(std::string(mesh)), Mesh_error, expected_msg);

  mesh = root_dir + mesh_dir + "missing_meshformat_end_spec.msh";
  error_msg = "[" + mesh +
              "] The following specifiers "
              "are missing in the mesh file:\n$EndMeshFormat\n";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(std::string(mesh)), Mesh_error, expected_msg);
}

BOOST_AUTO_TEST_CASE(meshformat_content_existence) {
  std::string mesh(root_dir + mesh_dir + "missing_meshformat_content.msh");
  error_msg = "[" + mesh +
              "] Missing content between "
              "$MeshFormat specifiers.";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(std::string(mesh)), Mesh_error, expected_msg);
}

BOOST_AUTO_TEST_CASE(physicalnames_content_existence) {
  BOOST_CHECK_THROW(
      Check_mesh(
          root_dir + mesh_dir + "missing_physicalnames_content.msh"),
      Mesh_error);
}

BOOST_AUTO_TEST_CASE(multiple_specifiers) {
  std::string mesh(root_dir + mesh_dir + "missing_multiple_specs.msh");
  error_msg = "[" + mesh +
              "] The following specifiers are missing "
              "in the mesh "
              "file:\n$PhysicalNames\n$Entities\n$Nodes\n$Elements\n$"
              "EndPhysicalNames\n$EndEntities\n$EndNodes\n$"
              "EndElements\n";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(std::string(mesh)), Mesh_error, expected_msg);
}

BOOST_AUTO_TEST_SUITE_END();
} // namespace Mesh
