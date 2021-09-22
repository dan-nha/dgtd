#include "../../../src/spatial_solver/mesh/check_mesh.h"
#include "../../../src/tools/custom_errors.h"

#include <boost/test/unit_test.hpp>
#include <string>

namespace Mesh {

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
  error_msg =
      "[" + root_dir + mesh_dir + "empty.msh] File seems to be empty.";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(root_dir + mesh_dir + "empty.msh"),
      Mesh_error,
      expected_msg);
}

BOOST_AUTO_TEST_CASE(empty_lines) {
  error_msg =
      "[" + root_dir + mesh_dir +
      "empty_lines.msh] Detected empty line in mesh file. (line 1)";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(root_dir + mesh_dir + "empty_lines.msh"),
      Mesh_error,
      expected_msg);
}

BOOST_AUTO_TEST_CASE(meshformat_specifiers) {
  error_msg =
      "[" + root_dir + mesh_dir +
      "missing_meshformat_begin_spec.msh] The following specifiers are "
      "missing in the mesh file:\n$MeshFormat\n";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(
          root_dir + mesh_dir + "missing_meshformat_begin_spec.msh"),
      Mesh_error,
      expected_msg);

  error_msg = "[" + root_dir + mesh_dir +
              "missing_meshformat_end_spec.msh] The following specifiers "
              "are missing in the mesh file:\n$EndMeshFormat\n";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(root_dir + mesh_dir + "missing_meshformat_end_spec.msh"),
      Mesh_error,
      expected_msg);
}

BOOST_AUTO_TEST_CASE(meshformat_content) {
  error_msg = "[" + root_dir + mesh_dir +
              "missing_meshformat_content.msh] Missing content between "
              "$MeshFormat specifiers.";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(root_dir + mesh_dir + "missing_meshformat_content.msh"),
      Mesh_error,
      expected_msg);
}

BOOST_AUTO_TEST_CASE(physicalnames_content) {
  error_msg = "[" + root_dir + mesh_dir +
              "missing_physicalnames_content.msh] Missing content between "
              "$PhysicalNames specifiers.";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(
          root_dir + mesh_dir + "missing_physicalnames_content.msh"),
      Mesh_error,
      expected_msg);
}

BOOST_AUTO_TEST_CASE(multiple_specifiers) {
  error_msg =
      "[" + root_dir + mesh_dir +
      "missing_multiple_specs.msh] The following specifiers are missing "
      "in the mesh "
      "file:\n$PhysicalNames\n$Entities\n$Nodes\n$Elements\n$"
      "EndPhysicalNames\n$EndEntities\n$EndNodes\n$"
      "EndElements\n";
  BOOST_CHECK_EXCEPTION(
      Check_mesh(root_dir + mesh_dir + "missing_multiple_specs.msh"),
      Mesh_error,
      expected_msg);
}
} // namespace Mesh
