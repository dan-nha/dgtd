# make spatial_solver library
file(GLOB_RECURSE spatial_solver_sources ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")
add_library(spatial_solver SHARED ${spatial_solver_sources})

install(TARGETS spatial_solver DESTINATION lib/internal)
