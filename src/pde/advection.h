#ifndef ADVECTION_H
#define ADVECTION_H

#include "pde.h"

class Advection: public Pde{
  public:
    Advection(const double flux_prefactor);

    inline std::tuple<double, double> get_boundary_conditions(
      const arma::mat &fields,
      const double time) const override {
      return {-sin(flux_prefactor*time), 0.};
    };

    arma::mat get_initial_values(
      const arma::mat &phys_node_coords) const override;

    arma::mat get_surface_fields(
      const arma::mat &fields,
      const double time,
      const std::vector<double> &geometric_factors,
      const arma::mat &lift_matrix,
      const double upwind_param) const override;

    inline double get_volume_flux_prefactor() const override {
      return flux_prefactor;
    };

    std::vector<std::tuple<double, double>> 
      get_surface_flux_prefactors(const size_t num_elems) const override;

  inline std::list<std::string> get_field_names() const override {
    return {"Advection"};
  };

  private:
    const double flux_prefactor;

};

#endif
