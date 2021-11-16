#include "jacobi_basis.h"

#include <boost/log/trivial.hpp>
#include <cmath>

namespace DG {

std::vector<double> Jacobi_basis::get_gauss_lobatto_nodes(
    const double alpha,
    const double beta,
    const size_t polynomial_order) const {

  if (polynomial_order < 1) {
    throw std::invalid_argument(
        std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
        ": "
        "To calculate the Gauss-Lobatto quadrature nodes the polynomial "
        "order must be >= 1.");
  } else {
    std::vector<double> gl_nodes(polynomial_order + 1);
    gl_nodes.front() = -1.0;
    for (size_t i{1}; i < polynomial_order; ++i) {
      gl_nodes[i] = this->get_gauss_jacobi_nodes(
          alpha + 1.0, beta + 1.0, polynomial_order - 2)(i - 1);
    }
    gl_nodes.back() = 1.0;
    return gl_nodes;
  }
}
//-------------------------------------------------------------------------
double Jacobi_basis::get_jacobi_polynomial(
    const double alpha,
    const double beta,
    const size_t polynomial_order,
    const double position) const {
  if (alpha <= -1.0) {
    throw std::invalid_argument(
        std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
        ": "
        "Parameter alpha must greater than -1.");
  }
  if (beta <= -1.0) {
    throw std::invalid_argument(
        std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
        ": "
        "Parameter beta must be greater than -1.");
  }
  if ((alpha + beta) == -1.0) {
    throw std::invalid_argument(
        std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
        ": "
        "To calculate the Jacobi polynomials the "
        "sum alpha+beta must not be -1.");
  };

  arma::vec jacobi_poly(polynomial_order + 1);

  if (polynomial_order == 0) {
    return get_recurrence_initiation(
               alpha, beta, polynomial_order, position)
        .front();
  } else if (polynomial_order == 1) {
    return get_recurrence_initiation(
               alpha, beta, polynomial_order, position)
        .back();
  } else {
    std::vector<double> initial_jacobi_poly{
        get_recurrence_initiation(alpha, beta, polynomial_order, position)};
    jacobi_poly[0] = initial_jacobi_poly.front();
    jacobi_poly[1] = initial_jacobi_poly.back();

    this->apply_reccurence(
        jacobi_poly, alpha, beta, polynomial_order, position);
  }
  return jacobi_poly[polynomial_order];
}
//-------------------------------------------------------------------------
double Jacobi_basis::get_jacobi_polynomial_gradient(
    const double alpha,
    const double beta,
    const size_t polynomial_order,
    const double position) const {

  if (polynomial_order == 0)
    return 0.0;

  return std::sqrt(
             double(polynomial_order) *
             (double(polynomial_order) + alpha + beta + 1.)) *
         this->get_jacobi_polynomial(
             alpha + 1., beta + 1., polynomial_order - 1, position);
}
//-------------------------------------------------------------------------
arma::vec Jacobi_basis::get_gauss_jacobi_nodes(
    const double alpha,
    const double beta,
    const size_t polynomial_order) const {
  arma::vec gj_nodes{arma::eig_sym(
      this->get_jacobi_polynomial_matrix(alpha, beta, polynomial_order))};

  if (gj_nodes.size() != polynomial_order + 1) {
    throw std::invalid_argument(
        std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
        ": "
        "The number of Gauss Jacobi quadrature nodes must be equal to "
        "polynomial_order+1.");
  };
  return gj_nodes;
}
//-------------------------------------------------------------------------
arma::mat Jacobi_basis::get_jacobi_polynomial_matrix(
    const double alpha,
    const double beta,
    const size_t polynomial_order) const {

  if (alpha <= -1.0) {
    throw std::invalid_argument(
        std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
        ": "
        "Parameter alpha must be >= 0.");
  }
  if (beta <= -1.0) {
    throw std::invalid_argument(
        std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
        ": "
        "Parameter beta must be >= 0.");
  }
  if (polynomial_order < 0) {
    throw std::invalid_argument(
        std::string{} + __FILE__ + ":" + std::to_string(__LINE__) +
        ": "
        "Polynomial order must be a positive integer.");
  }

  arma::vec diagonal{get_matrix_diagonal(alpha, beta, polynomial_order)};
  arma::vec subdiagonal{
      get_matrix_subdiagonal(alpha, beta, polynomial_order)};

  return get_tridiagonal_matrix(diagonal, subdiagonal, polynomial_order);
};
//----
arma::mat Jacobi_basis::get_tridiagonal_matrix(
    const arma::vec diagonal,
    const arma::vec subdiagonal,
    const size_t polynomial_order) const {

  arma::mat tridiag_matrix(
      polynomial_order + 1, polynomial_order + 1, arma::fill::zeros);
  tridiag_matrix.diag(0) = diagonal;

  if (polynomial_order > 0) {
    tridiag_matrix.diag(-1) = subdiagonal;
    tridiag_matrix.diag(1) = subdiagonal;
  }

  return tridiag_matrix;
}
//----
arma::vec Jacobi_basis::get_matrix_diagonal(
    const double alpha,
    const double beta,
    const size_t polynomial_order) const {

  arma::vec diagonal(polynomial_order + 1);
  diagonal(0) = (beta - alpha) / (2 + beta + alpha);

  for (size_t j{1}; j < polynomial_order + 1; ++j) {
    const double aux(2. * j + alpha + beta);
    diagonal(j) = (beta * beta - alpha * alpha) / aux / (aux + 2.);
  }

  return diagonal;
}
//----
arma::vec Jacobi_basis::get_matrix_subdiagonal(
    const double alpha,
    const double beta,
    const size_t polynomial_order) const {

  arma::vec subdiagonal(polynomial_order);

  for (size_t j{1}; j < polynomial_order + 1; ++j) {
    const double aux(2. * j + alpha + beta);
    subdiagonal(j - 1) = 2. / aux *
                         std::sqrt(
                             j * (j + alpha + beta) * (j + alpha) *
                             (j + beta) / (aux + 1.) / (aux - 1.));
  };

  return subdiagonal;
}
//-------------------------------------------------------------------------
std::vector<double> Jacobi_basis::get_recurrence_initiation(
    const double alpha,
    const double beta,
    const size_t polynomial_order,
    const double position) const {

  std::vector<double> jacobi_ini;

  const double aux0{
      std::pow(2., alpha + beta + 1.) / (alpha + beta + 1.) *
      std::tgamma(alpha + 1.) * std::tgamma(beta + 1.) /
      std::tgamma(alpha + beta + 1.)};
  jacobi_ini.push_back(1. / std::sqrt(aux0));

  if (polynomial_order >= 1) {
    const double aux1{
        (alpha + 1.) * (beta + 1.) / (alpha + beta + 3.) * aux0};
    jacobi_ini.push_back(
        ((alpha + beta + 2.) * position + (alpha - beta)) / 2. /
        std::sqrt(aux1));
  }

  return jacobi_ini;
}
//----
void Jacobi_basis::apply_reccurence(
    arma::vec &jacobi_polynomial,
    const double alpha,
    const double beta,
    const size_t polynomial_order,
    const double position) const {

  double a_old{
      2. / (2. + alpha + beta) *
      sqrt((alpha + 1) * (beta + 1) / (alpha + beta + 3))};
  double a_new, b_new, aux;

  for (size_t i{0}; i < polynomial_order - 1; ++i) {
    aux = 2. * double(i + 1) + alpha + beta;

    a_new = 2. / (aux + 2.) *
            std::sqrt(
                (double(i) + 2.) * (double(i) + 2. + alpha + beta) *
                (double(i) + 2. + alpha) * (double(i) + 2. + beta) /
                (aux + 1.) / (aux + 3.));

    b_new = -(alpha * alpha - beta * beta) / (aux * aux + 2 * aux);

    jacobi_polynomial[i + 2] =
        1. / a_new *
        (-a_old * jacobi_polynomial[i] +
         (position - b_new) * jacobi_polynomial[i + 1]);

    a_old = a_new;
  }
}
} // namespace DG
