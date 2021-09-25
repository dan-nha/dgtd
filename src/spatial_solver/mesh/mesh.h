#ifndef MESH_H
#define MESH_H

#include <list>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>

namespace DG::Mesh {

/// Struct for different entities occuring in the mesh
const struct Entity_type {
  size_t point, curve, surface, volume;
} constexpr Entity = {0, 1, 2, 3};

/// Struct for the physical groups
const struct Physical_group_type {
  size_t contour, region;
} constexpr Physical_group = {0, 1};

/**
 * Lookup table for physical groups and their respective entities depending
 * on the physical dimension
 */
inline static std::map<size_t, size_t>
physgroup_entity_table(size_t dimension) {
  switch (dimension) {
  case 1:
    return {
        {Physical_group.contour, Entity.point},
        {Physical_group.region, Entity.curve}};
  case 2:
    return {
        {Physical_group.contour, Entity.curve},
        {Physical_group.region, Entity.surface}};
  case 3:
    return {
        {Physical_group.contour, Entity.surface},
        {Physical_group.region, Entity.volume}};
  default:
    throw std::invalid_argument("Unknown physical dimension.");
  }
};

/**
 * Lookup table for entities and their respective physical groups depending
 * on the physical dimension
 */
inline static std::map<size_t, size_t>
entity_physgroup_table(size_t dimension) {
  switch (dimension) {
  case 1:
    return {
        {Entity.point, Physical_group.contour},
        {Entity.curve, Physical_group.region}};
  case 2:
    return {
        {Entity.curve, Physical_group.contour},
        {Entity.surface, Physical_group.region}};
  case 3:
    return {
        {Entity.surface, Physical_group.contour},
        {Entity.volume, Physical_group.region}};
  default:
    throw std::invalid_argument("Unknown physical dimension.");
  }
};

/// List of specifiers for each section in the mesh file
const std::list<std::string> specifier_list(
    {"$MeshFormat", "$PhysicalNames", "$Entities", "$Nodes", "$Elements"});

/**
 * @brief Structure for the different element geometries. Initial values as
 * described in gmsh (see https://gmsh.info/doc/texinfo/gmsh.html).
 * Gmsh has even more element types, which you might want to add (see
 * https://gitlab.onelab.info/gmsh/gmsh/blob/master/Common/GmshDefines.h
 * for more information Make sure that the added element types in the
 * Element_type enum class have the same unsigned ineger value as in gmsh.
 * One might want to think about replacing the enum class by a container in
 * the long term.
 */
const struct Element_type {
  size_t point;
  size_t line_2nodes;
  size_t triangle_3nodes;
  size_t tetrahedron_4nodes;
  size_t max_elem_type;
} constexpr Element{15, 1, 2, 4, 140};

static std::map<size_t, size_t> dimension_to_finite_elem_type{
    {1, Element.line_2nodes},
    {2, Element.triangle_3nodes},
    {3, Element.tetrahedron_4nodes}};

} // namespace DG::Mesh

#endif
