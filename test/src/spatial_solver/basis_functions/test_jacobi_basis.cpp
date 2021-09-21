#include "../../../../src/spatial_solver/basis_functions/jacobi_basis.h"

#include <boost/test/unit_test.hpp>
#include <cmath>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

namespace DG_solver {
Jacobi_basis jacobi;

BOOST_AUTO_TEST_CASE(gauss_jacobi_nodes, *utf::tolerance(1e-16)) {
  // hard-coded (tabulated) values were taken from
  // https://en.wikipedia.org/wiki/Gaussian_quadrature

  const double alpha(0), beta(0);
  const size_t polynomial_order(2);

  arma::vec nodes(
      jacobi.get_gauss_lobatto_nodes(alpha, beta, polynomial_order));
  BOOST_TEST(nodes[0] == -1, tt::tolerance(1e-15));
  BOOST_TEST(nodes[1] == 0, tt::tolerance(1e-15));
  BOOST_TEST(nodes[2] == 1, tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_CASE(jacobi_polynomials, *utf::tolerance(1e-16)) {

  /**
   * Reference values taken from Matlab nodal DGTD code by Hesthaven
   * and Warburton \cite hesthaven2008nodal (e.g. Jacobi.m)
   * (see also https://github.com/tcew/nodal-dg/tree/master/Codes1.1)
   **/

  double polynomial_order(0);
  BOOST_TEST(
      jacobi.get_jacobi_polynomial(1, 2, polynomial_order, 1.) ==
          0.866025403784439,
      tt::tolerance(1e-15));
  BOOST_TEST(
      jacobi.get_jacobi_polynomial(1, 2, polynomial_order, 8.23522) ==
          0.866025403784439,
      tt::tolerance(1e-15));

  polynomial_order = 1;
  BOOST_TEST(
      jacobi.get_jacobi_polynomial(1, 3, polynomial_order, 1.) ==
          1.479019945774904,
      tt::tolerance(1e-15));
  BOOST_TEST(
      jacobi.get_jacobi_polynomial(1, 3, polynomial_order, 5.379) ==
          11.193962459597360,
      tt::tolerance(1e-15));

  polynomial_order = 2;
  BOOST_TEST(
      jacobi.get_jacobi_polynomial(1, 3, polynomial_order, 1.) == 2.25,
      tt::tolerance(1e-15));
  BOOST_TEST(
      jacobi.get_jacobi_polynomial(3, 1, polynomial_order, 0.14) ==
          0.095400000000000,
      tt::tolerance(1e-14));

  polynomial_order = 3;
  BOOST_TEST(
      jacobi.get_jacobi_polynomial(4, 2, polynomial_order, 0.2527) ==
          -0.720585645170871,
      tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_CASE(jacobi_polynomial_gradient, *utf::tolerance(1e-16)) {
  BOOST_TEST(
      jacobi.get_jacobi_polynomial_gradient(0, 0, 1, 2.) ==
          1.224744871391589,
      tt::tolerance(1e-15));
  BOOST_TEST(
      jacobi.get_jacobi_polynomial_gradient(1, 3, 1, 0.2657) ==
          2.218529918662356,
      tt::tolerance(1e-15));
  BOOST_TEST(
      jacobi.get_jacobi_polynomial_gradient(1, 3, 2, 2.) == 18.375,
      tt::tolerance(1e-14));
}
} // namespace DG_solver
