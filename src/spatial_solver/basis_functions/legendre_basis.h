#ifndef LEGENDRE_BASIS_H
#define LEGENDRE_BASIS_H

#include "jacobi_basis.h"
#include "polynomial.h"

#include <cstddef>
#include <vector>

namespace DG {
/**
 * @brief Create Legendre polynomials, their gradients, and
 * quadrature nodes for the spatial solver.
 *
 * Within this class I create the Legendre polynomials, their gradients,
 * and the corresponding Gaussian quadrature nodes from the Jacobi
 * polynomials and Gauss-Lobatto nodes, respectively.
 */
class Legendre_basis : public Jacobi_basis, public Polynomial {

public:
  Legendre_basis(){};

  /**
   * @brief The Gaussian quadrature nodes for the Legendre polynomials are
   * a special case of the Gauss-Lobatto nodes of the Jacobi_basis class,
   * where the parameters of the Gauss-Lobatto nodes are set to
   * \f$\alpha=\beta=0\f$.
   */
  std::vector<double> get_quad_nodes(const size_t polynomial_order) const;

  /**
   * @brief The Legendre polynomials are derived from the Jacobi
   * polynomials by setting the Jacobi polynomial parameters to
   * \f$\alpha=\beta=0\f$.
   */
  double get_polynomial(
      const size_t polynomial_order,
      const double position) const;

  /**
   * @brief The Legendre polynomial gradient is derived from the Jacobi
   * polynomial gradient by setting the Jacobi polynomial gradient
   * parameters to \f$\alpha=\beta=0\f$.
   */
  double get_polynomial_gradient(
      const size_t polynomial_order,
      const double position) const;
};
} // namespace DG
#endif
