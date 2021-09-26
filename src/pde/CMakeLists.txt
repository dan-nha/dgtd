# make pde library
file(GLOB_RECURSE pde_sources ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")
add_library(pde SHARED ${pde_sources})

install(TARGETS pde DESTINATION lib/internal)
