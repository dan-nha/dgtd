# make executable

file(GLOB_RECURSE test_sources ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")

add_definitions(-DDGTD_ROOT=\"${PROJECT_SOURCE_DIR}\")

set(test_exec test_${PROJECT_NAME})
add_executable(${test_exec} ${test_sources})

# link libraries to executable
target_link_libraries(${test_exec} PRIVATE
  #spatial_solver
  temporal_solver
  #pde
  #tools
  ${Boost_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
)

install(TARGETS ${test_exec} DESTINATION bin)
