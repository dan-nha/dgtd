#ifndef PDE_H
#define PDE_H

#include <armadillo>
#include <list>
#include <map>
#include <string>
#include <tuple>

class Pde {
public:
  Pde(){};

  /**
   * @brief Calculate the boundary condition of a field at a certain time
   */
  virtual std::tuple<double, double> get_boundary_conditions(
      const arma::mat &fields,
      const double time) const = 0;
  
  /**
   * @brief Calculate initial values from the physical DG scheme quadrature
   * node coordinates
   */
  virtual arma::mat get_initial_values(
      const arma::mat &phys_node_coords) const = 0;
  
  /**
   * @brief Create the DG scheme for the PDE, where I distinguish between
   * volume fields and surface fields. The latter connect the individual
   * solution per element to a global solution via the numerical flux.
   */
  arma::mat get_spatial_scheme(
      const arma::mat &fields,
      const double time,
      const std::vector<double> &geometric_factors,
      const arma::mat &diff_matrix,
      const arma::mat &lift_matrix,
      const double upwind_param=0) const;

  /**
   * @brief Calculate the volume field, i.e.
   * \f$ \underline{\mathbf v}_h \equiv \boldsymbol{\mathcal{D}} \cdot
   * \underline{\mathbf f}_h \f$, where \f$ \boldsymbol{\mathcal{D}} \f$ is
   * the element independent \f$(N \times N)\f$-differentiation matrix.
   *
   * @param[in] fluxes Flux of the pervious time step's field
   * @param[in] geometric_factors Geometric factors of each element to
   * convert physical integrals from [a, b] to reference integrals from
   * [-1,1].
   * @param[in] diff_matrix Differentiation matrix for reference integrals
   * (element length independent)
   *
   * @return Matrix of volume fields, where each column represents an
   * element
   */
  arma::mat get_volume_fields(
      const arma::mat &fields,
      const std::vector<double> &geometric_factors,
      const arma::mat &diff_matrix) const;
  
  /**
   * @brief Calculate the surface field, i.e. \f$ \underline{\mathbf s}_h =
   * (\underline{s}_h^1, \dots \underline{s}_h^K) \f$, where each column is
   * given by <br> \f$ \underline{s}^k_h = \boldsymbol{\mathcal{L}}(x)
   * \left[f_h^k(u_h^k(x)) - f^*(u_h^k(x))
   * \right]_{x_\mathrm{L}^k}^{x_\mathrm{R}^k} \f$ <br> with the left
   * boundary \f$ x_\mathrm{L}^k \f$ and the right boundary \f$
   * x_\mathrm{R}^k \f$ of an element \f$ k\f$, the \f$(N \times 2)\f$-lift
   * matrix \f$\mathcal L\f$, and the numerical flux \f$f^*\f$. The lift
   * matrix at the left boundary shall be given by its first column,
   * whereas the lift matrix at the right boundary shall be given by its
   * second one. The numerical flux is given as in @cite hesthaven2008nodal
   * (p. 20 (chapter 2.1), p. 25 (chapter 2.2)).
   *
   * @param[in] fields Fields of the previous time step
   * @param[in] time Current time which is passed to the boundary
   * condition contained in the volume field calculation
   * @param[in] geometric_factors Geometric factors of each element to
   * convert physical integrals from [a, b] to reference integrals from
   * [-1,1].
   * @param[in] lift_matrix Lift matrix for reference integrals (element
   * length independent)
   * @param[in] upwind_param Upwind parameter \f$ \in [0,1]\f$, corresponds
   * to \f$(1-\alpha)\f$ in Hesthaven and Warburton's textbook
   * @cite hesthaven2008nodal (p.25, chapter 2.2).
   *
   * @return Matrix of surface fields, where each column represents an
   * element
   */
  virtual arma::mat get_surface_fields(
      const arma::mat &fields,
      const double time,
      const std::vector<double> &geometric_factors,
      const arma::mat &lift_matrix,
      const double upwind_param) const;
 
  inline arma::mat
  get_flux(const arma::mat &fields, const double flux_prefactor) const {
    return flux_prefactor * fields;
  }

  /// @brief Calculate field jumps at every interface within a region
  std::vector<std::tuple<double, double>> get_field_jumps(
      const arma::mat &fields,
      const std::tuple<double, double> boundary_conditions) const;
  
  /**
   * @brief Calculate left and right flux jump of an element multiplied by
   * the lift matrix and the corresponding prefactors depending on the
   * element's left or right boundary
   *
   * @param[in] prefactors Vector of tuples, where each tuple consists of
   * the prefactors at the left and right boundary of each element
   * @param[in] lift_matrix Lift matrix for (element length independent)
   * reference integrals
   * @param[in] geometri_factors Vector of geometric factors (Jacobians)
   * @param[in] upwind_param Upwind parameter \f$ \in [0,1]\f$, corresponds
   * to \f$(1-\alpha)\f$ in Hesthaven and Warburton's textbook
   *
   * @return Lifted flux jump matrix
   */
  arma::mat get_lifted_jumps(
      const std::vector<std::tuple<double, double>> &field_jumps,
      const std::vector<std::tuple<double, double>> &flux_prefactors,
      const arma::mat &lift_matrix,
      const std::vector<double> &geometric_factors,
      const double upwind_param) const;
 
  virtual double get_volume_flux_prefactor() const = 0;
  virtual std::vector<std::tuple<double, double>>
    get_surface_flux_prefactors(const size_t num_elems) const = 0;
  virtual std::list<std::string> get_field_names() const = 0;
};
#endif
