#define BOOST_TEST_MODULE test_temporal_solver

#include "../../../src/temporal_solver/low_storage_runge_kutta.h"
#include <boost/test/unit_test.hpp>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

namespace TD_solver {
Low_storage_runge_kutta lsrk_solver;
/**
 * Reference values taken from Matlab nodal DGTD code by Hesthaven
 * and Warburton \cite hesthaven2008nodal .
 * (see also https://github.com/tcew/nodal-dg/tree/master/Codes1.1)
 **/

//-------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(get_butcher_coeffs, *utf::tolerance(1e-16)) {

  auto [a, b, c] = lsrk_solver.get_butcher_coeffs(4, 5);
  BOOST_REQUIRE(a[0] == 0.);
  BOOST_TEST(a[1] == -0.417890474499852, tt::tolerance(1e-15));
  BOOST_TEST(a[2] == -1.192151694642677, tt::tolerance(1e-15));
  BOOST_TEST(a[3] == -1.697784692471528, tt::tolerance(1e-15));
  BOOST_TEST(a[4] == -1.514183444257156, tt::tolerance(1e-15));

  BOOST_TEST(b[0] == 0.149659021999229, tt::tolerance(1e-15));
  BOOST_TEST(b[1] == 0.379210312999627, tt::tolerance(1e-15));
  BOOST_TEST(b[2] == 0.822955029386982, tt::tolerance(1e-15));
  BOOST_TEST(b[3] == 0.699450455949122, tt::tolerance(1e-15));
  BOOST_TEST(b[4] == 0.153057247968152, tt::tolerance(1e-15));

  BOOST_REQUIRE(c[0] == 0.);
  BOOST_TEST(c[1] == 0.149659021999229, tt::tolerance(1e-15));
  BOOST_TEST(c[2] == 0.370400957364205, tt::tolerance(1e-15));
  BOOST_TEST(c[3] == 0.622255763134443, tt::tolerance(1e-15));
  BOOST_TEST(c[4] == 0.958282130674690, tt::tolerance(1e-15));

  BOOST_CHECK_THROW(lsrk_solver.get_butcher_coeffs(14, 4), std::exception);
}
//-------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(evolve_in_time, *utf::tolerance(1e-16)) {

  auto test_odes = [](arma::mat x, double t) {
    arma::mat odes(x.n_rows, x.n_cols);
    for (size_t row(0); row < x.n_rows; ++row) {
      for (size_t column(0); column < x.n_cols; ++column) {
        odes(row, column) =
            (5 * t * t - x(row, column)) / exp(x(row, column) + t);
      }
    }
    return odes;
  };
  const size_t num_stages(5);
  const auto [a, b, c] = lsrk_solver.get_butcher_coeffs(4, num_stages);
  const double time(7.);
  const double dt(0.1);
  const arma::mat ini_vals({0});

  const arma::mat first_evolution = lsrk_solver.evolve_in_time(
      test_odes, a, b, c, num_stages, time, dt, ini_vals);

  BOOST_TEST(
      first_evolution.front() == 0.021330398714844, tt::tolerance(1e-13));
}
} // namespace TD_solver
