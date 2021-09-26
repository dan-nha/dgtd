#include "../../src/pde/advection.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(advection);

Advection advection(2*M_PI);

BOOST_AUTO_TEST_CASE(boundary_conditions) {
  const arma::mat fields({{1,2,3}});
  const double time(0.25);
  const auto [left_bc, right_bc] = 
    advection.get_boundary_conditions(fields, time);
  BOOST_TEST(left_bc == -1);
  BOOST_TEST(right_bc == 0);
}

BOOST_AUTO_TEST_CASE(initial_values) {
  const arma::mat nodes({{M_PI/2., 3, 5}});
  const arma::mat ini_vals(advection.get_initial_values(nodes));
  BOOST_TEST(ini_vals(0,0) == 1);
  BOOST_TEST(ini_vals(0,1) == 0);
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

BOOST_AUTO_TEST_SUITE_END();
