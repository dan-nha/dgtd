#ifndef DGTD_SOLVER_H
#define DGTD_SOLVER_H

namespace DG::TD {

template <class Pde, class Basis, class TD_solver> class Dgtd_solver {
public:
  Dgtd_solver(
      const std::string &mesh_name,
      const size_t polynomial_order,
      const double end_time);

  bool is_field_name_valid(const std::string &field_name) const;
};
} // namespace DG::TD

#endif
