cmake_minimum_required(VERSION 2.8)

set(COMPILER intel CACHE STRING "Compiler suite")
include(${CMAKE_CURRENT_SOURCE_DIR}/${COMPILER}.cmake)
project (ljmc C Fortran)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_VERBOSE_MAKEFILE 1)

include_directories(/opt/local/include)
link_directories(/opt/local/lib)

add_executable(ljmc main vgwspb.f90 dlsode.f rs.f)

target_link_libraries(ljmc ${LINKFLAGS} -lgsl ${LAPACK})
