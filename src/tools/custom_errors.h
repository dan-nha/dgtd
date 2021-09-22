#include <exception>
#include <string>

/// @brief Custom mesh error class for errors in mesh file produced by gmsh
class Mesh_error: public std::exception {

  /// Error message
  std::string error_message;

  public:
  Mesh_error( const std::string& custom_error_msg, const std::string& filename)
  : error_message("[" + filename + "]\n" + custom_error_msg) 
  {};

  /**
   * @brief Custom throw message
   * @return The error message argument, i.e. the error message itself
   */
  virtual const char* what() const throw(){
    return error_message.c_str();
  }
};

/// @brief Custom error class for unimplemented methods
class Not_implemented: public std::exception {

  std::string error_message;

  public:
  Not_implemented( const std::string& custom_error_msg )
    : error_message( custom_error_msg ) 
  {};

  /**
   * @brief Custom throw message
   * @return The error message argument, i.e. the error message itself
   */
  virtual const char* what() const throw() {
    return error_message.c_str();
  }
};
