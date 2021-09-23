#include "../../../../src/spatial_solver/mesh/import_mesh_data.h"
#include "../../../../src/spatial_solver/mesh/mesh.h"
#include "../../../../src/spatial_solver/mesh/mesh_geometry.h"
#include "../../../src/tools/custom_errors.h"

#include <boost/test/unit_test.hpp>

namespace DG_solver::Mesh {

std::string error_msg;
bool expected_msg(const Mesh_error &ex) {
  BOOST_CHECK_EQUAL(ex.what(), error_msg);
  return true;
};

const std::string root_dir(DGTD_ROOT);
const std::string mesh_dir("/test/spatial_solver/mesh/test_meshes/");

//-------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(import_mesh_1d);

Import_mesh_data example(root_dir + mesh_dir + "example.msh");
Import_mesh_data line(root_dir + mesh_dir + "line.msh");

BOOST_AUTO_TEST_CASE(dimension) { BOOST_TEST(line.get_dimension() == 1); }

BOOST_AUTO_TEST_CASE(physical_groups) {
  BOOST_TEST(example.import_gmsh_physical_groups()[1] == 0);
  BOOST_TEST(example.import_gmsh_physical_groups()[2] == 0);
  BOOST_TEST(example.import_gmsh_physical_groups()[3] == 1);
  BOOST_TEST(example.import_gmsh_physical_groups()[4] == 1);
  BOOST_TEST(example.import_gmsh_physical_groups()[5] == 1);
  BOOST_TEST(example.import_gmsh_physical_groups()[8] == 1);
}

BOOST_AUTO_TEST_CASE(entities) {
  BOOST_TEST(example.import_gmsh_entities(Entity.point)[5][0] == 2);
  BOOST_TEST(
      example.import_gmsh_entities(Entity.point)[2].empty() == true);
  BOOST_TEST(example.import_gmsh_entities(Entity.curve)[6][0] == 8);
}

BOOST_AUTO_TEST_CASE(nodes) {
  Import_mesh_data parametric(
      root_dir + mesh_dir + "fail_meshes/parametric_flag.msh");
  BOOST_CHECK_THROW(
      parametric.import_gmsh_nodes(Entity.point, 1), Mesh_error);
}

BOOST_AUTO_TEST_CASE(node_coordinates) {
  BOOST_TEST(
      line.import_gmsh_nodes(Entity.point, 1)[1][0] == 0); // x coordinate
  BOOST_TEST(
      line.import_gmsh_nodes(Entity.point, 1)[1][1] == 0); // y cooridnate
  BOOST_TEST(
      line.import_gmsh_nodes(Entity.point, 1)[1][2] == 0); // z cooridnate

  BOOST_TEST(
      line.import_gmsh_nodes(Entity.point, 2)[2][0] == 10); // x coordinate
  BOOST_TEST(
      line.import_gmsh_nodes(Entity.point, 2)[2][1] == 0); // y cooridnate
  BOOST_TEST(
      line.import_gmsh_nodes(Entity.point, 2)[2][2] == 0); // z cooridnate

  BOOST_TEST(
      line.import_gmsh_nodes(Entity.curve, 1)[3][0] == 1.491565972698457);
  BOOST_TEST(
      line.import_gmsh_nodes(Entity.curve, 1)[4][0] == 2.82608965819235);
  BOOST_TEST(
      line.import_gmsh_nodes(Entity.curve, 1)[7][0] == 6.04422962956386);
}

BOOST_AUTO_TEST_CASE(gmsh_element_node_tags) {
  BOOST_TEST(line.import_gmsh_elements(Entity.point, 1)[1][0] == 1);
  BOOST_TEST(line.import_gmsh_elements(Entity.point, 2)[2][0] == 2);
  BOOST_TEST(line.import_gmsh_elements(Entity.curve, 1)[3][0] == 1);
  BOOST_TEST(line.import_gmsh_elements(Entity.curve, 1)[3][1] == 3);
  BOOST_TEST(line.import_gmsh_elements(Entity.curve, 1)[13][0] == 12);
  BOOST_TEST(line.import_gmsh_elements(Entity.curve, 1)[13][1] == 2);
}

BOOST_AUTO_TEST_CASE(gmsh_element_type) {
  BOOST_TEST(line.import_gmsh_element_type(Entity.point, 1) == 15);
  BOOST_TEST(line.import_gmsh_element_type(Entity.point, 2) == 15);
  BOOST_TEST(line.import_gmsh_element_type(Entity.curve, 1) == 1);
}

BOOST_AUTO_TEST_SUITE_END();
//-------------------------------------------------------------------------
/**
 * The methods to import Gmsh files do not depend on the overall mesh
 * dimension. Hence, I have also implemented tests for 2D and 3D. Note,
 * however that in principle I will only write 1D code within this project.
 */
BOOST_AUTO_TEST_SUITE(import_mesh_2d_3d);

Import_mesh_data cylinder(root_dir + mesh_dir + "cylinder.msh");
Import_mesh_data sphere(root_dir + mesh_dir + "sphere.msh");

BOOST_AUTO_TEST_CASE(dimension) {
  BOOST_TEST(cylinder.get_dimension() == 2);
  BOOST_TEST(sphere.get_dimension() == 3);
}

BOOST_AUTO_TEST_CASE(wrong_dimension) {
  BOOST_CHECK_THROW(
      Import_mesh_data imd(
          root_dir + mesh_dir + "fail_meshes/wrong_dimension.msh"),
      Mesh_error);
}

BOOST_AUTO_TEST_CASE(physical_groups_2d) {
  BOOST_TEST(cylinder.import_gmsh_physical_groups()[1] == 1);
  BOOST_TEST(cylinder.import_gmsh_physical_groups()[2] == 1);
  BOOST_TEST(cylinder.import_gmsh_physical_groups()[3] == 1);
  BOOST_TEST(cylinder.import_gmsh_physical_groups()[4] == 1);
  BOOST_TEST(cylinder.import_gmsh_physical_groups()[5] == 2);
  BOOST_TEST(cylinder.import_gmsh_physical_groups()[6] == 2);
  BOOST_TEST(cylinder.import_gmsh_physical_groups()[7] == 2);
  BOOST_TEST(cylinder.import_gmsh_physical_groups()[8] == 2);
  BOOST_TEST(cylinder.import_gmsh_physical_groups()[9] == 2);
}

BOOST_AUTO_TEST_CASE(physical_groups_3d) {
  BOOST_TEST(sphere.import_gmsh_physical_groups()[1] == 2);
  BOOST_TEST(sphere.import_gmsh_physical_groups()[4] == 2);
  BOOST_TEST(sphere.import_gmsh_physical_groups()[5] == 2);
  BOOST_TEST(sphere.import_gmsh_physical_groups()[6] == 3);
  BOOST_TEST(sphere.import_gmsh_physical_groups()[7] == 3);
  BOOST_TEST(sphere.import_gmsh_physical_groups()[8] == 3);
  BOOST_TEST(sphere.import_gmsh_physical_groups()[9] == 3);
  BOOST_TEST(sphere.import_gmsh_physical_groups()[10] == 3);
  BOOST_TEST(sphere.import_gmsh_physical_groups()[11] == 3);
}

BOOST_AUTO_TEST_CASE(entities_2d) {
  BOOST_TEST(
      cylinder.import_gmsh_entities(Entity.point)[13].empty() == true);
  BOOST_TEST(cylinder.import_gmsh_entities(Entity.curve)[29][0] == 4);
  BOOST_TEST(cylinder.import_gmsh_entities(Entity.surface)[13][0] == 7);
}

BOOST_AUTO_TEST_CASE(entities_3d) {
  BOOST_TEST(
      sphere.import_gmsh_entities(Entity.point)[74].empty() == true);
  BOOST_TEST(sphere.import_gmsh_entities(Entity.curve)[77][0] == 2);
  BOOST_TEST(sphere.import_gmsh_entities(Entity.surface)[115][0] == 1);
  BOOST_TEST(sphere.import_gmsh_entities(Entity.volume)[29][0] == 9);
}

BOOST_AUTO_TEST_CASE(nodes_2d) {
  BOOST_TEST(
      cylinder.import_gmsh_nodes(Entity.point, 9)[9][0] ==
      200); // x coordinate
  BOOST_TEST(
      cylinder.import_gmsh_nodes(Entity.point, 9)[9][1] ==
      200); // y cooridnate
  BOOST_TEST(
      cylinder.import_gmsh_nodes(Entity.point, 9)[9][2] ==
      0); // z cooridnate

  BOOST_TEST(
      cylinder.import_gmsh_nodes(Entity.curve, 4)[51][0] ==
      59.99999999999997); // x coordinate
  BOOST_TEST(
      cylinder.import_gmsh_nodes(Entity.curve, 4)[51][1] ==
      150); // y cooridnate
  BOOST_TEST(
      cylinder.import_gmsh_nodes(Entity.curve, 4)[51][2] ==
      0); // z cooridnate

  BOOST_TEST(
      cylinder.import_gmsh_nodes(Entity.surface, 5)[349][0] ==
      340.9616587406966); // x coordinate
  BOOST_TEST(
      cylinder.import_gmsh_nodes(Entity.surface, 5)[349][1] ==
      -299.9999999999999); // y cooridnate
  BOOST_TEST(
      cylinder.import_gmsh_nodes(Entity.surface, 5)[349][2] ==
      0); // z cooridnate
}

BOOST_AUTO_TEST_CASE(import_gmsh_nodes_3d) {
  BOOST_TEST(
      sphere.import_gmsh_nodes(Entity.point, 2)[2][0] ==
      3.061616997868383e-15); // x cooridnate
  BOOST_TEST(
      sphere.import_gmsh_nodes(Entity.curve, 4)[95][2] ==
      107.1428571428571); // z cooridnate
  BOOST_TEST(
      sphere.import_gmsh_nodes(Entity.surface, 108)[2702][1] ==
      -5.684341886080801e-14); // y cooridnate
  BOOST_TEST(
      sphere.import_gmsh_nodes(Entity.volume, 29)[5307][0] ==
      139.0555123899853); // x cooridnate
}

BOOST_AUTO_TEST_CASE(gmsh_element_node_tags) {
  BOOST_TEST(cylinder.import_gmsh_elements(Entity.curve, 1)[3][0] == 23);
  BOOST_TEST(cylinder.import_gmsh_elements(Entity.curve, 1)[3][1] == 24);
  BOOST_TEST(
      cylinder.import_gmsh_elements(Entity.surface, 1)[89][2] == 205);
  BOOST_TEST(cylinder.import_gmsh_element_type(Entity.curve, 3) == 1);
  BOOST_TEST(cylinder.import_gmsh_element_type(Entity.surface, 1) == 2);
}

BOOST_AUTO_TEST_CASE(import_gmsh_element_node_tags) {
  BOOST_TEST(sphere.import_gmsh_elements(Entity.surface, 1)[7][2] == 748);
  BOOST_TEST(
      sphere.import_gmsh_elements(Entity.volume, 27)[39876][3] == 5150);
  BOOST_TEST(sphere.import_gmsh_element_type(Entity.surface, 1) == 2);
  BOOST_TEST(
      sphere.import_gmsh_element_type(Entity.volume, 27) ==
      Element.tetrahedron_4nodes);
}

BOOST_AUTO_TEST_SUITE_END();
} // namespace DG_solver::Mesh
