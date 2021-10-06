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
 * Advec1D.m
 * AdvecRHS1D.m
 * AdvecDriver1D.m: N=3, MeshGen1D(-3,6,3)
 */

const std::string root_dir(DGTD_ROOT);
const std::string mesh(root_dir + "/test/examples/three_elems.msh");
const size_t polynomial_order(3);
const double end_time(0.1);
const double dt_factor(0.75 * 0.5);
const double upwind_param(1.);

Dgtd_solver<Advection, Legendre_basis, Low_storage_runge_kutta>
    dgtd(mesh, polynomial_order, end_time, dt_factor, upwind_param);

Advection advection(2 * M_PI);

BOOST_AUTO_TEST_CASE(phys_node_coords, *utf::tolerance(1e-16)) {
  const arma::mat coords(dgtd.get_phys_node_coords());
  BOOST_TEST(coords(0, 0) == -3, tt::tolerance(1e-15));
  BOOST_TEST(coords(1, 0) == -2.17082039324994, tt::tolerance(1e-14));
  BOOST_TEST(coords(2, 0) == -0.829179606750063, tt::tolerance(1e-14));
  BOOST_TEST(coords(3, 0) == 0, tt::tolerance(1e-15));

  BOOST_TEST(coords(0, 1) == 0, tt::tolerance(1e-15));
  BOOST_TEST(coords(1, 1) == 0.829179606750063, tt::tolerance(1e-14));
  BOOST_TEST(coords(2, 1) == 2.17082039324994, tt::tolerance(1e-14));
  BOOST_TEST(coords(3, 1) == 3, tt::tolerance(1e-15));

  BOOST_TEST(coords(0, 2) == 3, tt::tolerance(1e-15));
  BOOST_TEST(coords(1, 2) == 3.82917960675006, tt::tolerance(1e-14));
  BOOST_TEST(coords(2, 2) == 5.17082039324994, tt::tolerance(1e-14));
  BOOST_TEST(coords(3, 2) == 6, tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_CASE(initial_values, *utf::tolerance(1e-16)) {
  const arma::mat coords(dgtd.get_phys_node_coords());
  const arma::mat ini_vals(advection.get_initial_values(coords));
  BOOST_TEST(ini_vals(0, 0) == -0.141120008059867, tt::tolerance(1e-14));
  BOOST_TEST(ini_vals(1, 0) == -0.825322025727965, tt::tolerance(1e-14));
  BOOST_TEST(ini_vals(2, 0) == -0.737377459323434, tt::tolerance(1e-14));
  BOOST_TEST(ini_vals(3, 0) == 0, tt::tolerance(1e-14));

  BOOST_TEST(ini_vals(0, 1) == 0, tt::tolerance(1e-14));
  BOOST_TEST(ini_vals(1, 1) == 0.737377459323434, tt::tolerance(1e-14));
  BOOST_TEST(ini_vals(2, 1) == 0.825322025727965, tt::tolerance(1e-14));
  BOOST_TEST(ini_vals(3, 1) == 0.141120008059867, tt::tolerance(1e-14));

  BOOST_TEST(ini_vals(0, 2) == 0.141120008059867, tt::tolerance(1e-14));
  BOOST_TEST(ini_vals(1, 2) == -0.634674278057035, tt::tolerance(1e-14));
  BOOST_TEST(ini_vals(2, 2) == -0.896747766176097, tt::tolerance(1e-14));
  BOOST_TEST(ini_vals(3, 2) == -0.279415498198926, tt::tolerance(1e-14));
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

BOOST_AUTO_TEST_CASE(initial_volume_fields, *utf::tolerance(5e-15)) {
  // Tested against -a*rx.*(Dr*u) in AdvecRHS1D.m
  const arma::mat coords(dgtd.get_phys_node_coords());
  const arma::mat ini_vals(advection.get_initial_values(coords));
  const std::vector<double> geo_factors(dgtd.get_geometric_factors());
  const arma::mat diff_matrix(
      Elementwise_operations<Legendre_basis>(polynomial_order)
          .get_diff_matrix());
  const arma::mat ini_volume_fields(
      advection.get_volume_fields(ini_vals, geo_factors, diff_matrix));
  BOOST_TEST(ini_volume_fields(0, 0) == 7.438566190284299);
  BOOST_TEST(ini_volume_fields(1, 0) == 2.975065528375445);
  BOOST_TEST(ini_volume_fields(2, 0) == -3.682489439168077);
  BOOST_TEST(ini_volume_fields(3, 0) == -7.448179281084407);

  BOOST_TEST(ini_volume_fields(0, 1) == -7.448179281084415);
  BOOST_TEST(ini_volume_fields(1, 1) == -3.682489439168079);
  BOOST_TEST(ini_volume_fields(2, 1) == 2.975065528375442);
  BOOST_TEST(ini_volume_fields(3, 1) == 7.438566190284294);

  BOOST_TEST(ini_volume_fields(0, 2) == 7.308717012932643);
  BOOST_TEST(ini_volume_fields(1, 2) == 4.316208298798116);
  BOOST_TEST(ini_volume_fields(2, 2) == -2.208095660804580);
  BOOST_TEST(ini_volume_fields(3, 2) == -7.280070146609996);
}

BOOST_AUTO_TEST_CASE(initial_surface_fields, *utf::tolerance(5e-15)) {
  // Tested against LIFT*(Fscale.*(du)) in AdvecRHS1D.m
  const arma::mat coords(dgtd.get_phys_node_coords());
  const arma::mat ini_vals(advection.get_initial_values(coords));
  const double time(0.);
  const std::vector<double> geo_factors(dgtd.get_geometric_factors());
  const arma::mat lift_matrix(
      Elementwise_operations<Legendre_basis>(polynomial_order)
          .get_lift_matrix());

  double upwind_param(0.);
  const arma::mat central_surface_fields(advection.get_surface_fields(
        ini_vals, time, geo_factors, lift_matrix, upwind_param));
  BOOST_TEST(central_surface_fields(0,0) == 2.364488429842194);
  BOOST_TEST(central_surface_fields(1,0) == -0.264357843056945);
  BOOST_TEST(central_surface_fields(2,0) == 0.264357843056944);
  BOOST_TEST(central_surface_fields(3,0) == -0.591122107460549);

  upwind_param = 1.;
   const arma::mat upwind_surface_fields(advection.get_surface_fields(
        ini_vals, time, geo_factors, lift_matrix, upwind_param)); 
  BOOST_TEST(upwind_surface_fields(0,0) == 4.728976859684387);
  BOOST_TEST(upwind_surface_fields(1,0) == -0.528715686113889);
  BOOST_TEST(upwind_surface_fields(2,0) == 0.528715686113889);
  BOOST_TEST(upwind_surface_fields(3,0) == -1.182244214921098);
}

/*
BOOST_AUTO_TEST_CASE(evolved_volume_fields, *utf::tolerance(1e-16)) {
  // Tested against -a*rx.*(Dr*u) in AdvecRHS1D.m
  const arma::mat coords(dgtd.get_phys_node_coords());
  const arma::mat ini_vals(advection.get_initial_values(coords));
  const std::vector<double> geo_factors(dgtd.get_geometric_factors());
  const arma::mat diff_matrix(
      Elementwise_operations<Legendre_basis>(polynomial_order)
          .get_diff_matrix());
  const arma::mat ini__volume_fields(
      advection.get_volume_fields(ini_vals, geo_factors, diff_matrix));
  const arma::mat evolved_volume_fields;
}
*/

/*
BOOST_AUTO_TEST_CASE(initial_spatial_scheme, *utf::tolerance(1e-15)) {
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
  const arma::mat ini_spatial_scheme(advection.get_spatial_scheme(
      ini_vals,
      time,
      geo_factors,
      diff_matrix,
      lift_matrix,
      upwind_param));

  BOOST_TEST(
      ini_spatial_scheme(0, 0) == 12.167543049968685,
      tt::tolerance(1e-14));
  BOOST_TEST(
      ini_spatial_scheme(1, 0) == 2.446349842261556, tt::tolerance(1e-14));
  BOOST_TEST(
      ini_spatial_scheme(2, 0) == -3.153773753054189,
      tt::tolerance(1e-14));
  BOOST_TEST(
      ini_spatial_scheme(3, 0) == -8.630423496005506,
      tt::tolerance(1e-14));

  BOOST_TEST(
      ini_spatial_scheme(0, 1) == -7.448179281084415,
      tt::tolerance(1e-14));
  BOOST_TEST(
      ini_spatial_scheme(1, 1) == -3.682489439168079,
      tt::tolerance(1e-14));
  BOOST_TEST(
      ini_spatial_scheme(2, 1) == 2.975065528375442, tt::tolerance(1e-14));
  BOOST_TEST(
      ini_spatial_scheme(3, 1) == 7.438566190284294, tt::tolerance(1e-14));

  BOOST_TEST(
      ini_spatial_scheme(0, 2) == 7.308717012932643, tt::tolerance(1e-14));
  BOOST_TEST(
      ini_spatial_scheme(1, 2) == 4.316208298798116, tt::tolerance(1e-14));
  BOOST_TEST(
      ini_spatial_scheme(2, 2) == -2.208095660804580,
      tt::tolerance(1e-14));
  BOOST_TEST(
      ini_spatial_scheme(3, 2) == -7.280070146609996,
      tt::tolerance(1e-14));
}
*/

/*
BOOST_AUTO_TEST_CASE(solution) {
  std::cout.precision(15);
  dgtd.get_solution(advection, 4, 5).raw_print(std::cout, "solution"); }
  */

BOOST_AUTO_TEST_SUITE_END();
} // namespace DGTD
