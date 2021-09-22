#ifndef MESH_H
#define MESH_H

#include <list>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>

namespace Spatial_solver::Mesh {

/// Struct for different entities occuring in the mesh
const struct Entity_type {
  /// @cond Doxygen_Suppress
  size_t point, curve, surface, volume;
  /// @endcond
} constexpr Entity = {0, 1, 2, 3};

/// Struct for the physical groups
const struct Physical_group_type {
  /// @cond Doxygen_Suppress
  size_t contour, region;
  /// @endcond
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
} // namespace Spatial_solver::Mesh

#endif
