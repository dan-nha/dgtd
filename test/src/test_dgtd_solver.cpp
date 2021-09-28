#include "../../src/dgtd_solver.h"
#include "../../src/pde/advection.h"
#include "../../src/spatial_solver/basis_functions/legendre_basis.h"
#include "../../src/temporal_solver/low_storage_runge_kutta.h"

#include <boost/test/unit_test.hpp>

namespace DGTD {

using namespace DG;
using namespace TD;

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(dgtd_solver);
/**
 * Advec1D.m
 * AdvecRHS1D.m
 * AdvecDriver1D.m: MeshGen1D(-3,6,3)
 */

const std::string root_dir(DGTD_ROOT);
const std::string mesh(root_dir+"/test/examples/three_elems.msh");
const size_t polynomial_order(3);
const double end_time(0.1);
const double dt_factor(0.75*0.5);
const double upwind_param(1.);

Dgtd_solver<Advection, Legendre_basis, Low_storage_runge_kutta>
  dgtd(mesh, polynomial_order, end_time, dt_factor, upwind_param);

BOOST_AUTO_TEST_CASE(phys_node_coords, *utf::tolerance(1e-16)) {
  const arma::mat coords(dgtd.get_phys_node_coords());
  BOOST_TEST(coords.front() == -3);
  BOOST_TEST(coords[1] == -2.170820393249937, tt::tolerance(1e-15));
  BOOST_TEST(coords.back() == 6);
}

BOOST_AUTO_TEST_CASE(min_node_dist, *utf::tolerance(1e-16)) {
 // Tested against xmin in Advec1D.m
 BOOST_TEST(dgtd.get_min_node_dist() == 0.829179606750063, tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_CASE(time_step, *utf::tolerance(1e-16)) {
  // Tested against dt in Advec1D.m
  BOOST_TEST(dgtd.get_time_step()==0.033333333333333, tt::tolerance(1e-14));
}

BOOST_AUTO_TEST_CASE(solution) {
  Advection advection(2*M_PI);
  std::cout << dgtd.get_solution(advection,4,5);
}

BOOST_AUTO_TEST_SUITE_END();
} // namespace DGTD
