#include "../../../src/spatial_solver/geometric_operations.h"

#include <boost/test/unit_test.hpp>

namespace DG {

const std::string root_dir(DGTD_ROOT);
const std::string mesh_dir("/test/src/spatial_solver/mesh/test_meshes/");

BOOST_AUTO_TEST_SUITE(geometric_operations);

Mesh::Process_mesh_data processed_mesh(root_dir + mesh_dir + "line.msh");
Geometric_operations line(processed_mesh);

BOOST_AUTO_TEST_CASE(geometric_factor) {
  BOOST_TEST(line.get_geometric_factor(3) == 2. / 1.491565972698457);
}

BOOST_AUTO_TEST_CASE(coord_conversion) {
  const double ref_coord(line.convert_phys_to_ref_coord(1., 3));
  BOOST_TEST(line.convert_ref_to_phys_coord(ref_coord, 3) == 1.);
}

BOOST_AUTO_TEST_CASE(min_dist) {
  const std::vector<double> coords({0.3,1,1.7,5,7});
  BOOST_TEST(line.get_min_distance(coords) == 0.7);
}

/**
 * get_min_node_dist(...) methods are tested in test_dgtd_solver.cpp
 */

BOOST_AUTO_TEST_SUITE_END();
} // namespace DG
