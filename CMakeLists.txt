cmake_minimum_required(VERSION 2.8)

# Change the default behavior by setting the environment variables CC, FC and
# BLA_VENDOR
# CC: icc, gcc, pgcc, etc.
# FC: ifort, gfortran, pgf90 etc.
# BLA_VENDOR: Intel10_lp64, ACML_MP, ACML_FM4_MP, ACML, etc.
set(CMAKE_C_COMPILER_INIT "icc")
set(CMAKE_Fortran_COMPILER_INIT "ifort")
set(BLA_VENDOR "mkl_parallel")

PROJECT(SCP C Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    set(CMAKE_C_FLAGS "-openmp")
    set(CMAKE_Fortran_FLAGS "-heap-arrays -openmp")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -ipo")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost -ipo")

    set(CMAKE_Fortran_FLAGS_DEBUG "-fpe0 -traceback -ftrapuv -check all  -warn unused -O0 -g")

    # cmake doesn't set these on Darwin??
    set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-O2 -g")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

    set(CMAKE_AR xiar)

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    set(CMAKE_C_FLAGS "-fopenmp")
    set(CMAKE_Fortran_FLAGS "-fopenmp -fimplicit-none")

    #set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    #set(CMAKE_C_FLAGS_RELEASE "-O3")
        set(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -ffpe-trap=invalid,zero,overflow,denormal -fimplicit-none -O0 -g")

elseif (CMAKE_Fortran_COMPILER_ID STREQUAL PGI OR CMAKE_C_COMPILER_ID STREQUAL PGI)
    set(CMAKE_C_FLAGS "-mp")
    set(CMAKE_Fortran_FLAGS "-mp")

        set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
        set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Ktrap=fp -Mbounds")
        set(CMAKE_C_FLAGS_RELEASE "-fast")
	set(BLA_VENDOR "ACML")
endif()

include(FindBLAS2.cmake)

include(FindCHOLMOD.cmake)

#set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

ADD_EXECUTABLE(scp  utils.f90 dlsode.f integrator.f90 lsode.f90 ek.f90
	vgw.f90 vgwfm.f90 cholmod_logdet.c csr_sort_ja.c reallocate.f90
	sparse.f90 ffvgw.f90 integration.f90 xyz.f90 scp.f90)

    
include_directories(scp ${CHOLMOD_INC})

target_link_libraries(scp ${CHOLMOD_LIBRARIES} ${BLAS_LIBRARIES})

if (APPLE)
  set_target_properties(scp PROPERTIES LINK_FLAGS -Wl,-no_pie)
endif()
