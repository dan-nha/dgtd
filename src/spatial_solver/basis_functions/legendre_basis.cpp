#include "legendre_basis.h"

#include <cmath>

namespace DG {
std::vector<double>
Legendre_basis::get_quad_nodes(const size_t polynomial_order) const {
  return this->get_gauss_lobatto_nodes(0, 0, polynomial_order);
};

double Legendre_basis::get_polynomial(
    const size_t polynomial_order,
    const double position) const {
  return this->get_jacobi_polynomial(0, 0, polynomial_order, position);
};

double Legendre_basis::get_polynomial_gradient(
    const size_t polynomial_order,
    const double position) const {
  return this->get_jacobi_polynomial_gradient(
      0, 0, polynomial_order, position);
};
} // namespace DG
