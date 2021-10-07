#include "../../src/dgtd_solver.h"
#include "../../src/pde/advection.h"
#include "../../src/spatial_solver/basis_functions/legendre_basis.h"
#include "../../src/spatial_solver/elementwise_operations.h"
#include "../../src/temporal_solver/low_storage_runge_kutta.h"
#include "../../src/tools/output.h"

#include <boost/test/unit_test.hpp>

#include <iomanip>

namespace DGTD {

using namespace DG;
using namespace TD;

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(dgtd_solver);
/**
 * Reference values taken from Matlab nodal DGTD code by Hesthaven
 * and Warburton \cite hesthaven2008nodal.
 * (see also https://github.com/tcew/nodal-dg/tree/master/Codes1.1)
 * More specifically, the code from AdvecDriver1D.m was used with the
 * following parameters N=3, MeshGen1D(0,9,3)
 */

const std::string root_dir(DGTD_ROOT);
const std::string
    mesh(root_dir + "/test/examples/advection/three_elems.msh");
const size_t polynomial_order(3);
const double end_time(0.1);
const double dt_factor(0.75 * 0.5);
const double upwind_param(1.);

Dgtd_solver<Advection, Legendre_basis, Low_storage_runge_kutta>
    dgtd(mesh, polynomial_order, end_time, dt_factor, upwind_param);

Advection advection(2 * M_PI);

BOOST_AUTO_TEST_CASE(phys_node_coords, *utf::tolerance(1e-15)) {
  const arma::mat coords(dgtd.get_phys_node_coords());
  BOOST_TEST(coords(0, 0) == 0);
  BOOST_TEST(coords(1, 0) == 0.829179606750063);
  BOOST_TEST(coords(2, 0) == 2.170820393249937);
  BOOST_TEST(coords(3, 0) == 3);

  BOOST_TEST(coords(0, 1) == 3);
  BOOST_TEST(coords(1, 1) == 3.829179606750063);
  BOOST_TEST(coords(2, 1) == 5.170820393249937);
  BOOST_TEST(coords(3, 1) == 6);

  BOOST_TEST(coords(0, 2) == 6);
  BOOST_TEST(coords(1, 2) == 6.829179606750063);
  BOOST_TEST(coords(2, 2) == 8.170820393249937);
  BOOST_TEST(coords(3, 2) == 9);
}

BOOST_AUTO_TEST_CASE(initial_values, *utf::tolerance(1e-14)) {
  const arma::mat coords(dgtd.get_phys_node_coords());
  const arma::mat ini_vals(advection.get_initial_values(coords));
  BOOST_TEST(ini_vals(0, 0) == 0);
  BOOST_TEST(ini_vals(1, 0) == 0.737377459323434);
  BOOST_TEST(ini_vals(2, 0) == 0.825322025727965);
  BOOST_TEST(ini_vals(3, 0) == 0.141120008059867);

  BOOST_TEST(ini_vals(0, 1) == 0.141120008059867);
  BOOST_TEST(ini_vals(1, 1) == -0.634674278057035);
  BOOST_TEST(ini_vals(2, 1) == -0.896747766176097);
  BOOST_TEST(ini_vals(3, 1) == -0.279415498198926);

  BOOST_TEST(ini_vals(0, 2) == -0.279415498198926);
  BOOST_TEST(ini_vals(1, 2) == 0.519268086800104);
  BOOST_TEST(ini_vals(2, 2) == 0.950225093987128);
  BOOST_TEST(ini_vals(3, 2) == 0.412118485241757);
}

BOOST_AUTO_TEST_CASE(geo_factors, *utf::tolerance(1e-16)) {
  // Tested against Fscale in StartUp1D.m, where the Jacobian is calculated
  // in GeometricFactors1D.m
  for (const auto geo_factor : dgtd.get_geometric_factors()) {
    BOOST_TEST(geo_factor == 0.666666666666666, tt::tolerance(1e-15));
  }
}

BOOST_AUTO_TEST_CASE(min_node_dist, *utf::tolerance(1e-16)) {
  // Tested against xmin in Advec1D.m
  BOOST_TEST(
      dgtd.get_min_node_dist() == 0.829179606750063, tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_CASE(time_step, *utf::tolerance(1e-16)) {
  // Tested against dt in Advec1D.m
  BOOST_TEST(
      dgtd.get_time_step() == 0.033333333333333, tt::tolerance(1e-14));
}

BOOST_AUTO_TEST_CASE(initial_spatial_scheme, *utf::tolerance(1e-14)) {
  // Tested against rhsu in AdvecRHS1D.m
  const arma::mat coords(dgtd.get_phys_node_coords());
  const arma::mat ini_vals(advection.get_initial_values(coords));
  const std::vector<double> geo_factors(dgtd.get_geometric_factors());
  const arma::mat lift_matrix(
      Elementwise_operations<Legendre_basis>(polynomial_order)
          .get_lift_matrix());
  const arma::mat diff_matrix(
      Elementwise_operations<Legendre_basis>(polynomial_order)
          .get_diff_matrix());
  const double time(0.);
  const arma::mat spatial_scheme(advection.get_spatial_scheme(
      ini_vals,
      time,
      geo_factors,
      diff_matrix,
      lift_matrix,
      upwind_param));

  BOOST_TEST(spatial_scheme(0, 0) == -7.448179281084415);
  BOOST_TEST(spatial_scheme(1, 0) == -3.682489439168079);
  BOOST_TEST(spatial_scheme(2, 0) == 2.975065528375442);
  BOOST_TEST(spatial_scheme(3, 0) == 7.438566190284294);

  BOOST_TEST(spatial_scheme(0, 1) == 7.308717012932643);
  BOOST_TEST(spatial_scheme(1, 1) == 4.316208298798116);
  BOOST_TEST(spatial_scheme(2, 1) == -2.208095660804580);
  BOOST_TEST(spatial_scheme(3, 1) == -7.280070146609996);

  BOOST_TEST(spatial_scheme(0, 2) == -7.022970724074257);
  BOOST_TEST(spatial_scheme(1, 2) == -4.863538219981340);
  BOOST_TEST(spatial_scheme(2, 2) == 1.396930743569629);
  BOOST_TEST(spatial_scheme(3, 2) == 6.975863449453324);
}

BOOST_AUTO_TEST_SUITE_END();
} // namespace DGTD
