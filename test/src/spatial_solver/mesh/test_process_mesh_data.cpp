#include "../../../../src/spatial_solver/mesh/mesh.h"
#include "../../../../src/spatial_solver/mesh/process_mesh_data.h"

#include <boost/test/unit_test.hpp>

namespace DG::Mesh {

const std::string root_dir(DGTD_ROOT);
const std::string mesh_dir("/test/src/spatial_solver/mesh/test_meshes/");

//-------------------------------------------------------------------------
BOOST_AUTO_TEST_SUITE(process_mesh);

Process_mesh_data line(root_dir + mesh_dir + "line.msh");
Process_mesh_data example(root_dir + mesh_dir + "example.msh");

BOOST_AUTO_TEST_CASE(finite_elems) {
  std::vector finite_elems(line.get_finite_elems(2));
  for (size_t i(0); i < 11; ++i) {
    BOOST_TEST(finite_elems[i] == i + 3);
  }

  const std::vector finite_elems_by_name(
      line.get_finite_elems("the_only_region"));
  BOOST_CHECK_EQUAL_COLLECTIONS(
      finite_elems_by_name.begin(),
      finite_elems_by_name.end(),
      finite_elems.begin(),
      finite_elems.end());

  BOOST_CHECK_THROW(
      line.get_finite_elems("invalid_region"), std::invalid_argument);

  BOOST_CHECK_THROW(
      line.get_finite_elems("outer_bc"), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(finite_elem_coords) {
  {
    const auto [left_coord, right_coord] = line.get_elem_coords(3);
    BOOST_TEST(left_coord == 0);
    BOOST_TEST(right_coord == 1.491565972698457);
  }

  {
    const auto [left_coord, right_coord] = line.get_elem_coords(5);
    BOOST_TEST(left_coord == 2.82608965819235);
    BOOST_TEST(right_coord == 4.020105179930608);
  }

  {
    const auto [left_coord, right_coord] = line.get_elem_coords(13);
    BOOST_TEST(left_coord == 9.509679649033021);
    BOOST_TEST(right_coord == 10);
  }
}

BOOST_AUTO_TEST_CASE(ordered_regions) {
  std::vector ordered_regions(example.get_ordered_regions());
  BOOST_TEST(ordered_regions[0] == 3);
  BOOST_TEST(ordered_regions[1] == 4);
  BOOST_TEST(ordered_regions[2] == 5);
  BOOST_TEST(ordered_regions[3] == 8);
}

BOOST_AUTO_TEST_CASE(ordered_finite_elems) {
  const std::vector<size_t> elem_tags({9, 5, 13, 3});
  const std::vector ordered_elem_tags(line.get_ordered_elems(elem_tags));
  BOOST_TEST(ordered_elem_tags[0] == 3);
  BOOST_TEST(ordered_elem_tags[1] == 5);
  BOOST_TEST(ordered_elem_tags[2] == 9);
  BOOST_TEST(ordered_elem_tags[3] == 13);
}

BOOST_AUTO_TEST_CASE(min_elem_size) {
  BOOST_TEST(line.get_min_elem_size() == 0.4903203509669787);
}

BOOST_AUTO_TEST_CASE(regional_min_elem_size) {
  BOOST_TEST(line.get_min_elem_size(2) == 0.4903203509669787);
  BOOST_TEST(
      line.get_min_elem_size("the_only_region") ==
      line.get_min_elem_size(2));

  BOOST_TEST(example.get_min_elem_size(3) == 1.0046848956386096);
  BOOST_TEST(
      example.get_min_elem_size("left_pml") ==
      example.get_min_elem_size(3));

  BOOST_CHECK_THROW(
      line.get_min_elem_size("invalid_region"), std::invalid_argument);

  BOOST_CHECK_THROW(
      line.get_min_elem_size("outer_bc"), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(finite_elem_size) {
  BOOST_TEST(line.get_elem_size(3) == 1.491565972698457);
  BOOST_TEST(line.get_elem_size(5) == 1.1940155217382582);
  BOOST_TEST(line.get_elem_size(13) == 0.4903203509669787);
}

BOOST_AUTO_TEST_SUITE_END();
} // namespace DG::Mesh
