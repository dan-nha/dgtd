#include "polynomial.h"

namespace DG {
double
Polynomial::get_min_node_distance(const std::vector<double> quad_nodes) const {

  std::vector<double> dist;
  for (size_t i(1); i < quad_nodes.size(); ++i) {
    dist.push_back(fabs(quad_nodes[i] - quad_nodes[i - 1]));
  }
  return *std::min_element(dist.begin(), dist.end());
}
} // namespace DG
