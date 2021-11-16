#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include "dgtd_solver.h"
#include "pde/advection.h"
#include "spatial_solver/basis_functions/legendre_basis.h"
#include "spatial_solver/mesh/check_mesh.h"
#include "spatial_solver/mesh/process_mesh_data.h"
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
  DG::Mesh::Process_mesh_data processed_mesh(argv[1]);

  Input input(argv[2]);

  if (input.pde_name == "advection") {
    DGTD::Dgtd_solver<
        Advection,
        DG::Legendre_basis,
        TD::Low_storage_runge_kutta>
        dgtd(processed_mesh, input);

    Advection advection(input.material_params.front(), input.upwind_param);

    dgtd.get_solution(advection);
  }
}
