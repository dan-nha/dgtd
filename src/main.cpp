#include "dgtd_solver.h"
#include "spatial_solver/basis_functions/legendre_basis.h"
#include "spatial_solver/mesh/check_mesh.h"
#include "temporal_solver/low_storage_runge_kutta.h"
#include "pde/advection.h"
#include "tools/output.h"



#include <iostream>

void stream_welcome_message() {
  std::cout << "\033[31m"
    << "---------------------------------------\n"
    << "Welcome to DGTD \n\n"
    << "This is a one-dimensional \n"
    << "discontinuous Galerkin time-domain \n"
    << "(DGTD) solver.\n"
    << "---------------------------------------\n\033[0m"
    << std::endl;
}

int main(int argc, char* argv[]) {
  stream_welcome_message();

  Mesh::Check_mesh check_mesh(argv[1]);
  const size_t polynomial_order(3);
  const double end_time(0.1);
  const double dt_factor(0.75 * 0.5);
  const double upwind_param(1.);

  DGTD::Dgtd_solver<Advection, DG::Legendre_basis, TD::Low_storage_runge_kutta>
    dgtd(argv[1], polynomial_order, end_time, dt_factor, upwind_param);
  Advection advection(2 * M_PI);
  dgtd.get_solution(advection, 4, 5);
}

