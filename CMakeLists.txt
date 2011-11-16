cmake_minimum_required(VERSION 2.8)

if(INTEL)
	message("Intel!")
endif()
if(GCC)
	message("GCC!")
endif()

PROJECT(OH Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
find_package(LAPACK REQUIRED)

ADD_EXECUTABLE(OH utils.f90 lsode.f90 vgw.f90 vgwfm.f90 dlsode.f rs.f xyz.f90 mainvars.f90 kubo.f90 correlations.f90 main.f90)

target_link_libraries(OH ${LAPACK_LIBRARIES})

set(LINK_FLAGS -Wl,-no_pie)

if (LAPACK_LINKER_FLAGS)
	set(LINK_FLAGS ${LINK_FLAGS} ${LAPACK_LINKER_FLAGS})
endif()
set_target_properties(OH PROPERTIES LINK_FLAGS ${LINK_FLAGS})
