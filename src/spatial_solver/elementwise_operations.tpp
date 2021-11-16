namespace DG {

template <class Basis>
Elementwise_operations<Basis>::Elementwise_operations(
    const size_t _polynomial_order)
    : polynomial_order(_polynomial_order),
      nodes(basis.get_quad_nodes(_polynomial_order)) {}
//-------------------------------------------------------------------------
template <class Basis>
arma::mat Elementwise_operations<Basis>::get_vandermonde_matrix() const {

  arma::mat vand_mat(this->nodes.size(), this->polynomial_order + 1);

  for (size_t node_idx{0}; node_idx < this->nodes.size(); ++node_idx) {
    for (size_t n{0}; n <= this->polynomial_order; ++n) {
      vand_mat(node_idx, n) =
          this->basis.get_polynomial(n, this->nodes[node_idx]);
    }
  }

  return vand_mat;
}
//-------------------------------------------------------------------------
template <class Basis>
arma::mat
Elementwise_operations<Basis>::get_diff_matrix() const {

  arma::mat diff_mat(
      this->polynomial_order + 1, this->polynomial_order + 1);
  if (this->polynomial_order == 0) {
    diff_mat.zeros();
  } else {
    diff_mat = this->get_grad_vandermonde_matrix() *
               arma::inv(this->get_vandermonde_matrix());
  }
  return diff_mat;
}
//-------------------------------------------------------------------------
template <class Basis>
arma::mat Elementwise_operations<Basis>::get_lift_matrix() const {

  arma::mat inverse_mass_matrix(
      this->get_vandermonde_matrix() *
      arma::trans(this->get_vandermonde_matrix()));

  arma::mat lift_matrix(this->polynomial_order + 1, 2, arma::fill::zeros);

  lift_matrix.col(0) = inverse_mass_matrix.col(0);
  lift_matrix.col(1) = inverse_mass_matrix.col(this->polynomial_order);

  return lift_matrix;
}
//-------------------------------------------------------------------------
template <class Basis>
arma::mat
Elementwise_operations<Basis>::get_grad_vandermonde_matrix() const {

  arma::mat grad_vand_mat(
      this->nodes.size(), this->polynomial_order + 1, arma::fill::zeros);

  for (size_t node_idx{0}; node_idx < this->nodes.size(); ++node_idx) {
    for (size_t n{0}; n <= this->polynomial_order; ++n) {
      grad_vand_mat(node_idx, n) =
          this->basis.get_polynomial_gradient(n, this->nodes[node_idx]);
    }
  }

  return grad_vand_mat;
}
} // namespace DG
