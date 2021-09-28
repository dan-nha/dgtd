#ifndef OUTPUT_H
#define OUTPUT_H

#include "../../lib/external/cnpy/cnpy.h"
#include <armadillo>
#include <string>

class Output {
  public:
    Output(){};
    void store_coords(const arma::mat phys_node_coords) const;
    void store_time(const double time) const;
    void store_fields(
        const std::string &field_name, 
        const arma::mat &fields) const;
};

#endif
