#ifndef MESH_GEOMETRY_H
#define MESH_GEOMETRY_H

#include<map>

namespace DG_solver::Mesh {
  /**
   * @brief Structure for the different element geometries. Initial values as
   * described in gmsh (see https://gmsh.info/doc/texinfo/gmsh.html).
   * Gmsh has even more element types, which you might want to add (see 
   * https://gitlab.onelab.info/gmsh/gmsh/blob/master/Common/GmshDefines.h for more
   * information Make sure that the added element types in the Element_type enum
   * class have the same unsigned ineger value as in gmsh. One might want to think
   * about replacing the enum class by a container in the long term.
   */
  const struct Element_type { 
    /// @cond Doxygen_Suppress
      size_t point
      , line_2nodes
      , triangle_3nodes
      , tetrahedron_4nodes
      , line_3nodes_2nd_order
      , triangle_6nodes_2nd_order
      , tetrahedron_10nodes_2nd_order
      , max_elem_type
      ;
    ///@endcond
  } constexpr Element{15, 1, 2, 4, 8, 9, 11, 140};


  /// Map finite element geometry to physical dimension
  static std::map<size_t, size_t> elem_type_to_dimension {
      { Element.line_2nodes, 1 }
    , { Element.triangle_3nodes, 2 }
    , { Element.triangle_6nodes_2nd_order, 2 }
    , { Element.tetrahedron_4nodes, 3 }
    , { Element.tetrahedron_10nodes_2nd_order, 3 }
  };
}
#endif
