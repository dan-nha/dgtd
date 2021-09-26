cmake_minimum_required(VERSION 3.7.2)
project(DGTD VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#----------------------------------------------------------------------------
# Compiler settings
#
# To build the oridnary MiniDGTD version type in 
# $ cmake .
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -std=c++17 -frtti")
#
# To build the debug MiniDGTD version type in
# $ cmake -DCMAKE_BUILD_TYPE=Debug .
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -ggdb -g3")

#----------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/internal")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
#----------------------------------------------------------------------------
# Libraries

# boost libs
#set(Boost_USE_STATIC_LIBS OFF) 
#set(Boost_USE_MULTITHREADED ON)  
#set(Boost_USE_STATIC_RUNTIME OFF)
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
find_package(Boost COMPONENTS log unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIRS}) 

# use armadillo without runtimelibrary, wrap blas and lapack
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
include_directories("${PROJECT_SOURCE_DIR}/lib/external/armadillo-9.800.1/include")

# npy libs
include_directories("${PROJECT_SOURCE_DIR}/lib/external/cnpy/include")
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_DONT_USE_WRAPPER -lblas -llapack -lz")

# Subdirectories
add_subdirectory(lib/external)
add_subdirectory(src)
add_subdirectory(test)
#add_subdirectory("doc")
