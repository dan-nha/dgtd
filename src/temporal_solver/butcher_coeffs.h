#include <vector>
#include <cmath>

namespace TD_solver {

  /**
   * @brief Butcher coefficients of the Runge-Kutta method. 
   * Each stage corresponds to a component of the three vectorial 
   * Butcher coefficients butcher_coeff1, butcher_coeff2, 
   * and butcher_coeff3.
   */
  struct Butcher_coeff { 
  	const std::vector<double> butcher_coeff1;
  	const std::vector<double> butcher_coeff2;
  	const std::vector<double> butcher_coeff3;
  }; 

  /**
   * @brief Low-storage five-stage fourth-order Runge-Kutta method.
   * Hesthaven and Warburton, Nodal Discontinuous Galerkin Methods, 
   * Chapter 3.4.
   */
  Butcher_coeff Lsrk_45_carpenter = {
    // butcher_coeff1 =
    { 0.0
    ,-567301805773.0  / 1357537059087.0
    ,-2404267990393.0 / 2016746695238.0
    ,-3550918686646.0 / 2091501179385.0
    ,-1275806237668.0 / 842570457699.0
    }
    // butcher_coeff2 =
    ,{ 1432997174477.0 / 9575080441755.0
     , 5161836677717.0 / 13612068292357.0
     , 1720146321549.0 / 2090206949498.0
     , 3134564353537.0 / 4481467310338.0
     , 2277821191437.0 / 14882151754819.0
    }
    // butcher_coeff3 =
    ,{ 0.0
     , 1432997174477.0 / 9575080441755.0
     , 2526269341429.0 / 6820363962896.0
     , 2006345519317.0 / 3224310063776.0
     , 2802321613138.0 / 2924317926251.0
    }
  };

}  
