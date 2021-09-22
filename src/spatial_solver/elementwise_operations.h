/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * CAUTION:
 * Due to instantiation issues the method implementation is stored in
 * a .tpp-file (not .cpp)
 * For more information, see
 * https://stackoverflow.com/questions/8752837/undefined-reference-to-template-class-constructor
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */

#ifndef ELEMENTWISE_OPERATIONS_H
#define ELEMENTWISE_OPERATIONS_H

#include <armadillo>
#include <string>

namespace DG_solver {
/**
 * @brief Here I compute the matrix operations needed for the DG (spatial)
 * scheme. Instead of computing the mass, stiffness and face matrix
 * directly, I aim at calculating the stiffness and face matrix times the
 * inverse mass matrix, resulting in the so-called "differentiation" and
 * "lift matrix", respectively @cite hesthaven2008nodal.<br> As the element
 * wise operations are functionally the same for any given basis or element
 * type, the basis and element type are chosen as template parameters.
 */
template <class Basis> class Elementwise_operations {

public:
  /**
   * @param[in] order corresponds to the order \f$n\f$ of the polynomial
   * \f$p^{\alpha,\beta}_n (x)\f$ with \f$n \geq 0\f$
   */
  Elementwise_operations(const size_t polynomial_order);

  /**
   * @brief Compute the Vandermonde matrix from some orthonormalized
   * polynomials. For example one can compute the Vandermonde matrix from
   * the orthonormalized Jacobi polynomials -- a function of two doubles
   * alpha and beta, an unsigned integer n, and a double pos.<br> The
   * function is adapted from [Vandermonde1D.m]
   * (https://github.com/tcew/nodal-dg/blob/master/Codes1.1/Codes1D/Vandermonde1D.m)
   * @cite hesthaven2008nodal
   *
   * @return Vandermonde matrix
   */
  arma::mat get_vandermonde_matrix() const;

  /**
   * @brief Compute the differentiation matrix from the Vandermonde
   * matrix.<br>
   * The function is adapted from [Dmatrix1D.m]
   * (https://github.com/tcew/nodal-dg/blob/master/Codes1.1/Codes1D/Dmatrix1D.m)
   * \cite hesthaven2008nodal (chapter 3.2, Dmatrix1D.m)
   */
  arma::mat get_differentiation_matrix() const;

  /**
   * @brief Compute the lift matrix
   * \f[
   *  \mathcal L = \mathcal M^{-1} \mathcal E
   * \f] where \f$\mathcal M = \mathcal V \mathcal V^{\mathrm{T}}\f$
   * denotes the mass matrix which is composed of the Vandermonde matrix
   * \f$ \mathcal V \f$ and its transposed. Furthermore, the matrix \f$
   * \mathcal E=(\hat{\mathbf e}_{1}, \hat{\mathbf e}_{N_\mathrm p})\f$
   * denotes a matrix which is composed of the unit vectors \f$\hat{\mathbf
   * e}_{1}
   * =(1,\underbrace{0,\dots,0}_{(N_\mathrm{p}-1)\text{-times}})^\mathrm{T}\f$
   * and
   * \f$\hat{\mathbf e}_{N_\mathrm{p}}
   * =(\underbrace{0,\dots,0}_{(N_\mathrm{p}-1)\text{-times}},1)^\mathrm{T}\f$
   * where each of the unit vectors has a length of
   * \f$N_\mathrm{p} = n + 1 \f$ (\f$n\f$: polynomial order).<br>
   * The function is adapted from [Lift1D.m]
   * (https://github.com/tcew/nodal-dg/blob/master/Codes1.1/Codes1D/Lift1D.m)
   * \cite hesthaven2008nodal (chapter 3.2, Lift1D.m)
   */
  arma::mat get_lift_matrix() const;

  /**
   * @brief Compute the gradient Vandermonde matrix from the gradient of
   * some orthonormalized polynomials at the nodal points generated in the
   * polynomial class. For example one can compute the gradient Vandermonde
   * matrix from the gradient of the orthonormalized Jacobi polynomials --
   * a function of two doubles alpha and beta, an unsigned integer n, and a
   * double pos. Within the gradient Vandermonde matrix routine the
   * position pos is set to the nodal points initalized in the
   * constructor.<br> The function is adapted from [GradVandermonde1D.m]
   * (https://github.com/tcew/nodal-dg/blob/master/Codes1.1/Codes1D/GradVandermonde1D.m)
   * \cite hesthaven2008nodal (chapter 3.2, GradVandermonde1D.m)
   *
   * @return Gradient Vandermonde matrix
   */
  arma::mat get_grad_vandermonde_matrix() const;

private:
  mutable Basis basis;
  const size_t polynomial_order;
  const arma::vec nodes;
};
} // namespace DG_solver

#include "elementwise_operations.tpp"

#endif
