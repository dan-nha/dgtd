#ifndef ADVECTION_H
#define ADVECTION_H

#include "pde.h"

class Advection : public Pde {
public:
  Advection(const double flux_prefactor, const double upwind_param);

  inline std::tuple<double, double> get_boundary_conditions(
      const arma::mat &fields,
      const double time) const override {
    return {-sin(advection_speed * time), 0.};
  };

  arma::mat
  get_initial_values(const arma::mat &phys_node_coords) const override;

  inline double get_volume_flux_prefactor() const override {
    return advection_speed;
  };

  std::vector<std::tuple<double, double>>
  get_surface_flux_prefactors(const size_t num_elems) const override;

  inline double get_upwind_param() const override {
    return this->upwind_param;
  };

  inline std::list<std::string> get_field_names() const override {
    return {"Advection"};
  };

private:
  const double advection_speed;
  const double upwind_param;
};

#endif
