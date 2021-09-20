# make temporal_solver library
file(GLOB_RECURSE temporal_solver_sources ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")
add_library(temporal_solver SHARED ${temporal_solver_sources})

install(TARGETS temporal_solver DESTINATION lib/internal)
