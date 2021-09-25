#include "../../src/spatial_solver/basis_functions/legendre_basis.h"
#include "../../src/spatial_solver/elementwise_operations.h"

#include <boost/test/unit_test.hpp>
#include <armadillo>
#include <cmath>

namespace utf = boost::unit_test;
namespace tt = boost::test_tools;

namespace DG {
BOOST_AUTO_TEST_CASE(vandermonde, *utf::tolerance(1e-16)) {
  /**
   * Reference values taken from Matlab nodal DGTD code by Hesthaven
   * and Warburton \cite hesthaven2008nodal (e.g. Vandermonde1D.m)
   * (see also https://github.com/tcew/nodal-dg/tree/master/Codes1.1)
   **/
  Elementwise_operations<Legendre_basis> eo(2);
  arma::mat vandermonde_matrix(eo.get_vandermonde_matrix());
  BOOST_TEST(
      vandermonde_matrix(0, 0) == 0.707106781186547, tt::tolerance(1e-15));
  BOOST_TEST(
      vandermonde_matrix(1, 0) == 0.707106781186547, tt::tolerance(1e-15));
  BOOST_TEST(
      vandermonde_matrix(2, 0) == 0.707106781186547, tt::tolerance(1e-15));
  BOOST_TEST(
      vandermonde_matrix(0, 1) == -1.224744871391589,
      tt::tolerance(1e-15));
  BOOST_TEST(vandermonde_matrix(1, 1) == 0.0, tt::tolerance(1e-15));
  BOOST_TEST(
      vandermonde_matrix(2, 1) == 1.224744871391589, tt::tolerance(1e-15));
  BOOST_TEST(
      vandermonde_matrix(0, 2) == 1.581138830084190, tt::tolerance(1e-15));
  BOOST_TEST(
      vandermonde_matrix(1, 2) == -0.790569415042095,
      tt::tolerance(1e-15));
  BOOST_TEST(
      vandermonde_matrix(2, 2) == 1.581138830084190, tt::tolerance(1e-15));
}

BOOST_AUTO_TEST_CASE(grad_vandermonde, *utf::tolerance(1e-16)) {
  /**
   * Reference values taken from Matlab nodal DGTD code by Hesthaven
   * and Warburton \cite hesthaven2008nodal (e.g. GradVandermonde1D.m)
   * (see also https://github.com/tcew/nodal-dg/tree/master/Codes1.1)
   **/
  {
    const double polynomial_order(1);
    Elementwise_operations<Legendre_basis> eo(polynomial_order);
    arma::mat grad_vandermonde_matrix(eo.get_grad_vandermonde_matrix());
    BOOST_TEST(grad_vandermonde_matrix(0, 0) == 0.0, tt::tolerance(1e-15));
    BOOST_TEST(grad_vandermonde_matrix(1, 0) == 0.0, tt::tolerance(1e-15));
    BOOST_TEST(
        grad_vandermonde_matrix(0, 1) == 1.224744871391589,
        tt::tolerance(1e-15));
    BOOST_TEST(
        grad_vandermonde_matrix(1, 1) == 1.224744871391589,
        tt::tolerance(1e-15));
  }

  {
    const double polynomial_order(3);
    Elementwise_operations<Legendre_basis> eo(polynomial_order);
    arma::mat grad_vandermonde_matrix(eo.get_grad_vandermonde_matrix());
    BOOST_TEST(grad_vandermonde_matrix(0, 0) == 0.0, tt::tolerance(1e-15));
    BOOST_TEST(grad_vandermonde_matrix(1, 0) == 0.0, tt::tolerance(1e-15));
    BOOST_TEST(grad_vandermonde_matrix(2, 0) == 0.0, tt::tolerance(1e-15));
    BOOST_TEST(grad_vandermonde_matrix(3, 0) == 0.0, tt::tolerance(1e-15));

    BOOST_TEST(
        grad_vandermonde_matrix(0, 1) == 1.224744871391589,
        tt::tolerance(1e-15));
    BOOST_TEST(
        grad_vandermonde_matrix(1, 1) == 1.224744871391589,
        tt::tolerance(1e-15));
    BOOST_TEST(
        grad_vandermonde_matrix(2, 1) == 1.224744871391589,
        tt::tolerance(1e-15));
    BOOST_TEST(
        grad_vandermonde_matrix(3, 1) == 1.224744871391589,
        tt::tolerance(1e-15));

    BOOST_TEST(
        grad_vandermonde_matrix(0, 2) == -4.743416490252569,
        tt::tolerance(1e-15));
    BOOST_TEST(
        grad_vandermonde_matrix(1, 2) == -2.121320343559642,
        tt::tolerance(1e-15));
    BOOST_TEST(
        grad_vandermonde_matrix(2, 2) == 2.121320343559642,
        tt::tolerance(1e-15));
    BOOST_TEST(
        grad_vandermonde_matrix(3, 2) == 4.743416490252569,
        tt::tolerance(1e-15));

    BOOST_TEST(
        grad_vandermonde_matrix(0, 3) == 11.224972160321824,
        tt::tolerance(1e-14));
    BOOST_TEST(grad_vandermonde_matrix(1, 3) == 0.0, tt::tolerance(1e-15));
    BOOST_TEST(grad_vandermonde_matrix(2, 3) == 0.0, tt::tolerance(1e-15));
    BOOST_TEST(
        grad_vandermonde_matrix(3, 3) == 11.224972160321824,
        tt::tolerance(1e-14));
  }
}
//---------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(differentiation_matrix, *utf::tolerance(1e-16)) {
  /**
   * Reference values taken from Matlab nodal DGTD code by Hesthaven
   * and Warburton \cite hesthaven2008nodal (e.g. Dmatrix1D.m)
   * (see also https://github.com/tcew/nodal-dg/tree/master/Codes1.1)
   **/
  {
    const double polynomial_order(1);
    Elementwise_operations<Legendre_basis> eo(polynomial_order);
    arma::mat differentiation_matrix(eo.get_differentiation_matrix());

    BOOST_TEST(differentiation_matrix(0, 0) == -0.5, tt::tolerance(1e-15));
    BOOST_TEST(differentiation_matrix(0, 1) == 0.5, tt::tolerance(1e-15));
    BOOST_TEST(differentiation_matrix(1, 0) == -0.5, tt::tolerance(1e-15));
    BOOST_TEST(differentiation_matrix(1, 1) == 0.5, tt::tolerance(1e-15));
  }

  {
    const double polynomial_order(2);
    Elementwise_operations<Legendre_basis> eo(polynomial_order);
    arma::mat differentiation_matrix(eo.get_differentiation_matrix());

    BOOST_TEST(differentiation_matrix(0, 0) == -1.5, tt::tolerance(1e-15));
    BOOST_TEST(differentiation_matrix(0, 1) == 2.0, tt::tolerance(1e-15));
    BOOST_TEST(differentiation_matrix(0, 2) == -0.5, tt::tolerance(1e-15));

    BOOST_TEST(differentiation_matrix(1, 0) == -0.5, tt::tolerance(1e-15));
    BOOST_TEST(differentiation_matrix(1, 1) == 0.0, tt::tolerance(1e-15));
    BOOST_TEST(differentiation_matrix(1, 2) == 0.5, tt::tolerance(1e-15));

    BOOST_TEST(differentiation_matrix(2, 0) == 0.5, tt::tolerance(1e-15));
    BOOST_TEST(differentiation_matrix(2, 1) == -2.0, tt::tolerance(1e-15));
    BOOST_TEST(differentiation_matrix(2, 2) == 1.5, tt::tolerance(1e-15));
  }
}
//---------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(lift_matrix, *utf::tolerance(1e-16)) {
  /**
   * Reference values taken from Matlab nodal DGTD code by Hesthaven and
   * Warburton \cite hesthaven2008nodal (e.g. Lift1D.m), where the
   * Gauss-Lobatto nodes are used as position points (see also
   * https://github.com/tcew/nodal-dg/tree/master/Codes1.1)
   **/

  {
    // Use Lift1D.m:
    // n=2; Np=n+1=3; Nfaces=2; Nfp=1; V=Vandermonde1D(n,JacobiGL(0,0,n));
    // Lift1D(Np, Nfaces, Nfp, V)
    const double polynomial_order(2);
    Elementwise_operations<Legendre_basis> eo(polynomial_order);
    arma::mat lift_matrix(eo.get_lift_matrix());

    BOOST_TEST(lift_matrix.n_cols == 2);
    BOOST_TEST(lift_matrix.n_rows == 3);

    BOOST_TEST(lift_matrix(0, 0) == 4.50, tt::tolerance(1e-15));
    BOOST_TEST(lift_matrix(1, 0) == -0.75, tt::tolerance(1e-15));
    BOOST_TEST(lift_matrix(2, 0) == 1.50, tt::tolerance(1e-15));

    BOOST_TEST(lift_matrix(0, 1) == 1.50, tt::tolerance(1e-15));
    BOOST_TEST(lift_matrix(1, 1) == -0.75, tt::tolerance(1e-15));
    BOOST_TEST(lift_matrix(2, 1) == 4.50, tt::tolerance(1e-15));
  }
  {
    // Use Lift1D.m:
    // n=3; Np=n+1=4; Nfaces=7; Nfp=1; V=Vandermonde1D(n,JacobiGL(0,0,n));
    // Lift1D(Np, Nfaces, Nfp, V)
    const double polynomial_order(3);
    Elementwise_operations<Legendre_basis> eo(polynomial_order);
    arma::mat lift_matrix(eo.get_lift_matrix());

    BOOST_TEST(lift_matrix.n_cols == 2);
    BOOST_TEST(lift_matrix.n_rows == 4);

    BOOST_TEST(
        lift_matrix(0, 0) == 8.000000000000004, tt::tolerance(1e-15));
    BOOST_TEST(
        lift_matrix(1, 0) == -0.894427190999917, tt::tolerance(1e-15));
    BOOST_TEST(
        lift_matrix(2, 0) == 0.894427190999917, tt::tolerance(1e-15));
    BOOST_TEST(
        lift_matrix(3, 0) == -2.000000000000003, tt::tolerance(1e-15));

    BOOST_TEST(
        lift_matrix(0, 1) == -2.000000000000003, tt::tolerance(1e-15));
    BOOST_TEST(
        lift_matrix(1, 1) == 0.894427190999917, tt::tolerance(1e-15));
    BOOST_TEST(
        lift_matrix(2, 1) == -0.894427190999917, tt::tolerance(1e-15));
    BOOST_TEST(
        lift_matrix(3, 1) == 8.000000000000004, tt::tolerance(1e-15));

    for (size_t i = 2; i < lift_matrix.n_cols; ++i) {
      BOOST_TEST(lift_matrix(0, i) == 0.0);
      BOOST_TEST(lift_matrix(1, i) == 0.0);
      BOOST_TEST(lift_matrix(2, i) == 0.0);
      BOOST_TEST(lift_matrix(3, i) == 0.0);
    }
  }
}
} // namespace DG
