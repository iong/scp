set(CMAKE_C_COMPILER icc)
set(CMAKE_Fortran_COMPILER ifort)
set(OPTFLAGS -O3 -mdynamic-no-pic -no-prec-div)
set(CMAKE_C_FLAGS_RELEASE ${OPTFLAGS})
set(CMAKE_C_FLAGS_DEBUG -O0 -g)
set(CMAKE_Fortran_FLAGS_RELEASE ${OPTFLAGS} -r8 -heap-arrays)
set(CMAKE_Fortran_FLAGS_DEBUG -O0 -g -r8 -heap-arrays -check bounds -ftrapuv -traceback)
set(LINKFLAGS -nofor-main)
enable_language(C)
enable_language(Fortran)
find_library(LAPACK NAMES mkl_sequential mkl_core PATHS $ENV{MKLROOT}/lib/em64t)


