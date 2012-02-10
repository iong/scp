cmake_minimum_required(VERSION 2.8)

if (NOT DEFINED ENV{FC})
    message(STATUS "NOT DEFINED ENV{FC}")
    set(ENV{CC} icc)
    set(ENV{FC} ifort)
endif()

PROJECT(SCP C Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

    set(MKL_LIBRARIES -mkl=parallel)
    set(LAPACK_LIBRARIES -mkl=parallel)
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -ipo -heap-arrays -openmp")
    set(CMAKE_Fortran_FLAGS_DEBUG "-fpe0 -traceback -ftrapuv -check all  -warn unused -O0 -g -heap-arrays -openmp")
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g -openmp")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost -ipo -openmp")
    set(CMAKE_AR xiar)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)

        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fimplicit-none")
        set(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -ffpe-trap=invalid,zero,overflow,denormal -fimplicit-none -O0 -g")
        set(CMAKE_C_FLAGS_RELEASE "-O3")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI OR CMAKE_C_COMPILER_ID STREQUAL PGI)

        set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
        set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Ktrap=fp -Mbounds")
        set(CMAKE_C_FLAGS_RELEASE "-fast")

endif()

include(FindMKL.cmake)

if (NOT DEFINED LAPACK_LIBRARIES)
	find_package(LAPACK REQUIRED)
endif()

#include(FindCHOLMOD.cmake)

#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

#Targets!

# Quantum Gibbs with VGW-SP
#ADD_EXECUTABLE(qgibbs utils.f90 dlsode.f vgw.f90 total_energy_vgw.f90 qgibbs.f90)

#target_link_libraries(qgibbs ${LAPACK_LIBRARIES})
#if (LAPACK_LINKER_FLAGS)
#    set_target_properties(qgibbs PROPERTIES LINK_FLAGS ${LAPACK_LINKER_FLAGS})
#endif()

add_subdirectory(nr)

include_directories(${nr_BINARY_DIR})
ADD_EXECUTABLE(clustergs  utils.f90 dlsode.f integrator.f90 lsode.f90 ek.f90
	vgw.f90 vgwfm.f90 integration.f90 clustergs.f90)

    
#include_directories(clustergs ${CHOLMOD_INC})

target_link_libraries(clustergs nr ${LAPACK_LIBRARIES})

add_executable(gaussmin gaussmin.f90 xyz.f90)
target_link_libraries(gaussmin nr)

add_executable(scpmin scpmin.f90 utils.f90 xyz.f90)
target_link_libraries(scpmin nr ${LAPACK_LIBRARIES})

if (APPLE)
  set_target_properties(clustergs PROPERTIES LINK_FLAGS -Wl,-no_pie)
endif()
