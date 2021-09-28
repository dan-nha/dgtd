#include "../../../src/pde/advection.h"

#include <boost/test/unit_test.hpp>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(advection);

Advection advection(2*M_PI);

BOOST_AUTO_TEST_CASE(boundary_conditions1) {
  const arma::mat fields({{1,2,3}});
  const double time(0.25);
  const auto [left_bc, right_bc] = 
    advection.get_boundary_conditions(fields, time);
  BOOST_TEST(left_bc == -1);
  BOOST_TEST(right_bc == 0);
}

BOOST_AUTO_TEST_CASE(boundary_conditions2, *utf::tolerance(1e-16)) {
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
  BOOST_TEST(left_bc == 2*M_PI);
  BOOST_TEST(left_bc == 2*M_PI);
}

BOOST_AUTO_TEST_CASE(field_name) {
  const std::string field_name(advection.get_field_names().front());
  BOOST_TEST(field_name == "Advection");
}

/**
 * get_surface_fields(...) is tested in test_dgtd_solver.cpp
 */

BOOST_AUTO_TEST_SUITE_END();
