#ifndef INPUT_H
#define INPUT_H

#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

class Input {  
  
  public:
    Input(const std::string& json_filename){
      pt::ptree root;
      pt::read_json(json_filename, root);

      end_time = root.get<double>("end_time");
      for (auto &&region_tree : root.get_child("regions")) {
        const pt::ptree &region_params = region_tree.second;

        pde_name = region_params.get<std::string>("pde");
        polynomial_order = region_params.get<size_t>("polynomial_order");
        runge_kutta_order = region_params.get<size_t>("runge_kutta_order");
        runge_kutta_stages =
        region_params.get<size_t>("runge_kutta_stages");
        dt_factor = region_params.get<double>("dt_factor");

        if (pde_name == "advection") {
          for (auto &&param_tree : region_params.get_child("parameters")) {
            const pt::ptree &pde_params = param_tree.second;
            material_params.push_back(
                pde_params.get<double>("advection_speed"));
            upwind_param = pde_params.get<double>("upwind_param");
          }
        }
      }
    };

    std::string pde_name;
    size_t polynomial_order;
    size_t runge_kutta_order;
    size_t runge_kutta_stages;
    double dt_factor;
    double end_time;
    double upwind_param;
    std::vector<double> material_params;
};

#endif
