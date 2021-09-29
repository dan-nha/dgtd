#include "output.h"

void Output::store_coords(const arma::mat phys_node_coords) const {
  cnpy::npy_save(
      "coordinates.npy",
      &phys_node_coords[0],
      {phys_node_coords.n_cols, phys_node_coords.n_rows},
      "w");
}
//-------------------------------------------------------------------------
void Output::store_time(const double time) const{
    cnpy::npy_save("times.npy", &time, {1}, "a");
  }
//-------------------------------------------------------------------------
void Output::store_fields(
    const std::string &field_name, 
    const arma::mat &fields) const {

  const std::string file_name(field_name + ".npy");
  cnpy::npy_save(
      file_name,
      &fields[0],
      {1,fields.n_cols, fields.n_rows},
      "a");
}
