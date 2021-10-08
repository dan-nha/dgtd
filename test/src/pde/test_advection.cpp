#include "../../../src/pde/advection.h"

#include <boost/test/unit_test.hpp>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(advection);

const double upwind_param(1.);
Advection advection(2*M_PI, upwind_param);

BOOST_AUTO_TEST_CASE(boundary_conditions1) {
  const arma::mat fields({{1,2,3}});
  const double time(0.25);
  const auto [left_bc, right_bc] = 
    advection.get_boundary_conditions(fields, time);
  BOOST_TEST(left_bc == -1);
  BOOST_TEST(right_bc == 0);
}

BOOST_AUTO_TEST_CASE(boundary_conditions2, *utf::tolerance(1e-16)) {
  /**
   * Tested against uin in AdvecRHS1D.m at the same points in time
   * Used AdvecDriver.m with parameters N=3 and MeshGen1D(-3,6,3)
   */
  std::vector<double> times(
      {0,
       0.004988634066641,
       0.012346698578807,
       0.020741858771148,
       0.031942737689156});
  arma::mat fields;

  std::vector<double> left_bc;
  for (const auto time: times) {
    left_bc.push_back(
      std::get<0>(advection.get_boundary_conditions(fields, time)));
  }
  BOOST_TEST(left_bc[0] == 0);
  BOOST_TEST(left_bc[1] == -0.031339379971237, tt::tolerance(1e-14));
  BOOST_TEST(left_bc[2] == -0.077498807531754, tt::tolerance(1e-13));
  BOOST_TEST(left_bc[3] == -0.129956336147177, tt::tolerance(1e-14));
  BOOST_TEST(left_bc[4] == -0.199357425830637, tt::tolerance(1e-13));
}

BOOST_AUTO_TEST_CASE(initial_values, *utf::tolerance(1e-16)) {
  const arma::mat nodes({{M_PI/2., M_PI, 0}});
  const arma::mat ini_vals(advection.get_initial_values(nodes));
  BOOST_TEST(ini_vals(0,0) == 1);
  BOOST_TEST(ini_vals(0,1) == 0, tt::tolerance(1e-15));
  BOOST_TEST(ini_vals(0,2) == 0);
}

BOOST_AUTO_TEST_CASE(volume_flux_prefactor) {
  BOOST_TEST(advection.get_volume_flux_prefactor() == 2*M_PI);
}

BOOST_AUTO_TEST_CASE(surface_flux_prefactor) {
  const size_t num_elems(1);
  const auto [left_bc, right_bc] = 
    advection.get_surface_flux_prefactors(num_elems).front();
  BOOST_TEST(left_bc == M_PI);
  BOOST_TEST(right_bc == M_PI);
}

BOOST_AUTO_TEST_CASE(field_jumps) {
  const double time(0.);
  const arma::mat field({{0, 2}, {5, 13}});
  const std::tuple<double, double> bc(
      advection.get_boundary_conditions(field, time));
  const std::vector<std::tuple<double, double>>
    jumps(advection.get_field_jumps(field, bc));

  BOOST_TEST(std::get<0>(jumps[0])==0);
  BOOST_TEST(std::get<1>(jumps[0])==3);
  BOOST_TEST(std::get<0>(jumps[1])==-3);
  BOOST_TEST(std::get<1>(jumps[1])==0);
}


BOOST_AUTO_TEST_CASE(lifted_field_jumps, *utf::tolerance(1e-15)) {
  const double time(0.);
  const arma::mat field({{0, 2}, {5, 13}});
  const std::tuple<double, double> bc(
      advection.get_boundary_conditions(field, time));
  const std::vector<std::tuple<double, double>>
    jumps(advection.get_field_jumps(field, bc));

  const arma::mat lift_mat({{2,3}});
  const std::vector<double> geo_factors({1,1});
  const std::vector<std::tuple<double, double>> prefactors(
      {{1,1}, {1,1}});

  double upwind_param(0.);
  const arma::mat central_jumps(advection.get_lifted_jumps(
        jumps, prefactors, lift_mat, geo_factors, upwind_param));
  BOOST_TEST(central_jumps.front() == 9);
  BOOST_TEST(central_jumps.back() == 6);

  upwind_param=0.1;
  const arma::mat upwind_jumps(advection.get_lifted_jumps(
          jumps, prefactors, lift_mat, geo_factors, upwind_param));
  BOOST_TEST(upwind_jumps.front() == -0.9);
  BOOST_TEST(upwind_jumps.back() == 0.6);
}

BOOST_AUTO_TEST_CASE(field_name) {
  const std::string field_name(advection.get_field_names().front());
  BOOST_TEST(field_name == "Advection");
}

/**
 * get_surface_fields(...) is tested in test_dgtd_solver.cpp by testing
 * get_spatial_scheme(...)
 */

BOOST_AUTO_TEST_SUITE_END();
