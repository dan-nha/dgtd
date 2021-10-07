#include "dgtd_solver.h"
#include "pde/advection.h"
#include "spatial_solver/basis_functions/legendre_basis.h"
#include "spatial_solver/mesh/check_mesh.h"
#include "temporal_solver/low_storage_runge_kutta.h"
#include "tools/input.h"

#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iostream>

namespace pt = boost::property_tree;

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

int main(int argc, char *argv[]) {
  stream_welcome_message();

  Mesh::Check_mesh check_mesh(argv[1]);

  Input input(argv[2]);

  if (input.pde_name == "advection") {
    DGTD::Dgtd_solver<
        Advection,
        DG::Legendre_basis,
        TD::Low_storage_runge_kutta>
        dgtd(
            argv[1],
            input.polynomial_order,
            input.end_time,
            input.dt_factor,
            input.upwind_param);

    Advection advection(input.advection_speed);
    
    dgtd.get_solution(
        advection, input.runge_kutta_order, input.runge_kutta_stages);
  }
}
