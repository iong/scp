cmake_minimum_required(VERSION 2.8)

add_executable(hsm-lj main.cpp)

target_link_libraries(hsm-lj ${BLAS_LIBRARIES})

if (APPLE)
  set_target_properties(hsm-lj PROPERTIES LINK_FLAGS -Wl,-no_pie)
endif()

