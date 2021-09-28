#ifndef OUTPUT_H
#define OUTPUT_H

#include "../../lib/external/cnpy/cnpy.h"
#include <armadillo>
#include <string>

class Output {
  void store_coords(const arma::mat phys_node_coords) const {
    cnpy::npy_save(
        "coordinates.npy",
        &phys_node_coords[0],
        {phys_node_coords.n_cols, phys_node_coords.n_rows},
        "w");
  }

  void store_time(const double time) const {
    cnpy::npy_save("times.npy", &time, {1}, "a");
  }

  void store_fields(const std::string &field_name, const arma::mat &fields)
      const {

    const std::string field_file_name(field_name + ".npy");
    cnpy::npy_save(
        field_file_name,
        &fields[0],
        {1, fields.n_cols, fields.n_rows},
        "a");
  }
};

#endif
