# add library target for cnpy 
add_library(cnpy STATIC "cnpy.cpp")
target_link_libraries(cnpy ${ZLIB_LIBRARIES})
