#ifndef JACOBI_BASIS_H
#define JACOBI_BASIS_H

#include <armadillo>
#include <cstddef>
#include <vector>

namespace DG_solver {
/**
 * @brief Create quadrature nodes for the spatial solver.
 *
 * Within this class we create the quadrature nodes for the spatial solver.
 * We use the Gauss Jacobi quadrature rule to create those nodes. A robust
 * and efficient way of implementing the Gauss-Jacobi quadrature rule is by
 * using the Golub & Welsch algorithm \cite golub1969calculation. This
 * algorithm takes the recurrence relation of the Jacobi polynomials and
 * creates a matrix vector problem of this relation. The resulting matrix
 * is of tridiagonal form. The matrix elements we have implemented here can
 * be found in \cite gil2007numerical [chapter 5.3, (5.109)]. One can show
 * that the eigenvalues of the tridiagonal matrix are identical to the
 * quadrature nodes. (The eigenvectors contain the quadrature weights.)
 */
class Jacobi_basis {

public:
  Jacobi_basis(){};

  /**
   * @brief Compute Gauss-Lobatto quadrature nodes.
   *
   * The outer nodes are the end points of the integration interval
   * \f$[-1,1]\f$.
   * The inner nodes are the Gauss-Jacobi quadrature nodes of
   * \f$P_{n-2}^{\alpha+1,\beta+1}\f$ computed with gauss_jacobi_nodes().
   * The function is adapted from [JacobiGL.m]
   * (https://github.com/tcew/nodal-dg/blob/master/Codes1.1/Codes1D/JacobiGL.m)
   * @cite hesthaven2008nodal.
   *
   * @note For polynomial order \f$n=1\f$ the only two nodes are the outer
   * nodes.
   */
  arma::vec get_gauss_lobatto_nodes(
      const double alpha,
      const double beta,
      const size_t polynomial_order) const;

  /**
   * @brief Compute the orthonormalized Jacobi polynomials
   *
   * We use the recurrence relation for orthonormalized Jacobi polynomials
   * as given in Hesthaven and Warburton \cite hesthaven2008nodal (Appendix
   * A). The Matlab reference code can be found under
   * /test/supplementary_material/hesthaven_warburton/JacobiP.m
   */
  double get_jacobi_polynomial(
      const double alpha,
      const double beta,
      const size_t polynomial_order,
      const double position) const;

  /**
   * @brief We calculate the gradient of the Jacobi polynomials by the
   * identity given in Hesthaven and Warburton \cite hesthaven2008nodal
   * (chap. 3.2). The Matlab reference code can be found under
   * /test/supplementary_material/hesthaven_warburton/GradJacobiP.m
   */
  double get_jacobi_polynomial_gradient(
      const double alpha,
      const double beta,
      const size_t order,
      const double position) const;

private:
  /**
   * @brief Compute Gauss-Jacobi quadrature nodes.
   *
   * Following the Golub-Welsch algorithm @cite golub1969calculation ,
   * the nodes are computed as the eigenvalues of the Jacobi matrix.
   * The form of the matrix elements can be found in @cite
   * gil2007numerical. The function is adapted from [JacobiGQ.m]
   * (https://github.com/tcew/nodal-dg/blob/master/Codes1.1/Codes1D/JacobiGQ.m)
   * @cite hesthaven2008nodal.
   *
   * @note The number of nodes is \f$n+1\f$.
   */
  arma::vec get_gauss_jacobi_nodes(
      const double alpha,
      const double beta,
      const size_t polynomial_order) const;

  /**
   * @brief Compute the symmetric tridiagonal Jacobi matrix for the
   * Gauss-Jacobi quadrature.
   *
   * Since the Jacobi matrix is symmetric and tridiagonal, all
   * nonzero elements can be stored in two vectors: the diagonal and the
   * subdiagonal (which is identical to the superdiagonal)
   */
  arma::mat get_jacobi_polynomial_matrix(
      const double alpha,
      const double beta,
      const size_t polynomial_order) const;

  arma::mat get_tridiagonal_matrix(
      const arma::vec diagonal,
      const arma::vec subdiagonal,
      const size_t polynomial_order) const;

  arma::vec get_matrix_diagonal(
      const double alpha,
      const double beta,
      const size_t polynomial_order) const;

  arma::vec get_matrix_subdiagonal(
      const double alpha,
      const double beta,
      const size_t polynomial_order) const;

  std::vector<double> get_recurrence_initiation(
      const double alpha,
      const double beta,
      const size_t polynomial_order,
      const double position) const;

  void apply_reccurence(
      arma::vec& jacobi_polynomial,
      const double alpha,
      const double beta,
      const size_t polynomial_order,
      const double position) const;
};
} // namespace DG_solver
#endif
