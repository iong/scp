cmake_minimum_required(VERSION 2.8)

set(COMPILER intel CACHE STRING "Compiler suite")
include(${CMAKE_CURRENT_SOURCE_DIR}/${COMPILER}.cmake)
project (ljmc C Fortran)

set(CMAKE_BUILD_TYPE Debug)

include_directories(/opt/local/include)
link_directories(/opt/local/lib)

add_executable(ljmc main vgwspb_H2_4G_Rc_Q_tau_SqrtMeff_Mar03.f dlsode.f)

target_link_libraries(ljmc ${LINKFLAGS} -lgsl ${LAPACK})
