#include "../../../../src/spatial_solver/basis_functions/legendre_basis.h"

#include <armadillo>
#include <boost/test/unit_test.hpp>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

namespace DG {
BOOST_AUTO_TEST_SUITE(basis_functions);
// hard-coded (tabulated) values were taken from
// https://en.wikipedia.org/wiki/Gaussian_quadrature
//
Legendre_basis gauss_legendre;

BOOST_AUTO_TEST_CASE(gauss_legendre_nodes, *utf::tolerance(1e-16)) {
  double polynomial_order(2);
  arma::vec nodes(gauss_legendre.get_quad_nodes(polynomial_order));
  BOOST_TEST(nodes[0] == -1.0);
  BOOST_TEST(nodes[1] == 0.0);
  BOOST_TEST(nodes[2] == 1.0);

  polynomial_order = 3;
  nodes = gauss_legendre.get_quad_nodes(polynomial_order);
  BOOST_TEST(nodes[0] == -1.0);
  BOOST_TEST(nodes[1] == -std::sqrt(1. / 5.));
  BOOST_TEST(nodes[2] == std::sqrt(1. / 5.));
  BOOST_TEST(nodes[3] == 1.0);
}

BOOST_AUTO_TEST_CASE(minimal_node_distance, *utf::tolerance(1e-16)) {
  arma::vec nodes({1, 4, 7, 12, 19});
  BOOST_TEST(gauss_legendre.get_min_node_distance(nodes) == 3);
}

BOOST_AUTO_TEST_SUITE_END();
} // namespace DG
