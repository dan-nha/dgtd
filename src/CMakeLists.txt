# Source files ana executable
file(GLOB_RECURSE dgtd_sources ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")
set(dgtd_exec ${PROJECT_NAME})
add_executable(${dgtd_exec} ${dgtd_sources})

# Link libraries to executable
target_link_libraries(${dgtd_exec} PRIVATE
    temporal_solver
    #spatial_solver
    #pde
    cnpy
    ${Boost_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${ZLIB}
    )
  

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

# Add subdirectories
#add_subdirectory("spatial_solver")
add_subdirectory("temporal_solver")
#add_subdirectory("pde")
add_subdirectory("tools")
