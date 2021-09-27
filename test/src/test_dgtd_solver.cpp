#include "../../src/dgtd_solver.h"
#include "../../src/pde/advection.h"
#include "../../src/spatial_solver/basis_functions/legendre_basis.h"
#include "../../src/temporal_solver/low_storage_runge_kutta.h"

#include <boost/test/unit_test.hpp>

namespace DGTD {
  using namespace DG;
  using namespace TD;

BOOST_AUTO_TEST_SUITE(dgtd_solver);

const std::string root_dir(DGTD_ROOT);
const std::string mesh_dir("/test/src/spatial_solver/mesh/test_meshes/");
const std::string mesh(root_dir+mesh_dir+"line.msh");
const size_t polynomial_order(3);
const double end_time(1.);

Dgtd_solver<Advection, Legendre_basis, Low_storage_runge_kutta>
  dgtd(mesh, polynomial_order, end_time);

BOOST_AUTO_TEST_CASE(phys_node_coords) {
  const arma::mat coords(dgtd.get_phys_node_coords());
  BOOST_TEST(coords.front() == 0);
  BOOST_TEST(coords.back() == 10);
}

BOOST_AUTO_TEST_SUITE_END();
} // namespace DGTD
