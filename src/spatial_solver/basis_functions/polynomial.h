#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <cstddef>
#include <vector>
#include <cmath>
#include <algorithm>

namespace DG {
/**
 * @brief General polynomial class
 *
 * Each class derived from the general polynomial class must have a
 * polynomial characterized by the order of the polynomial and the position
 * at which the polynomial is evaluated.<br>
 * Furthermore, any derived polynomial class must provide quadrature nodes
 * as well as the gradient of the polynomial.<br>
 */
class Polynomial {
public:
  Polynomial(){};

  virtual std::vector<double>
  get_quad_nodes(const size_t polynomial_order) const = 0;

  virtual double get_polynomial(
      const size_t polynomial_order,
      const double position) const = 0;

  virtual double get_polynomial_gradient(
      const size_t polynomial_order,
      const double position) const = 0;

  /**
   * @brief Depending on the type of finite element, the number of
   *quadrature nodes might vary. However, in 1D I always have a fixed
   *number of 1 as face quadrature node.
   **/
  inline size_t get_number_of_face_nodes() const { return 1; };

  /**
   * @brief Calculate the smallest distance between a given set of ordered
   * quadrature nodes. In 1D that order is from left to right.
   **/
  double get_min_node_distance(const std::vector<double> quad_nodes) const;
};
} // namespace DG
#endif
