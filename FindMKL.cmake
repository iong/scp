if (DEFINED MKL_LIBRARIES)
    return()
endif()

find_library(MKL_CORE NAMES mkl_core PATH_SUFFIXES intel64 HINTS $ENV{MKLROOT}/lib)
if (NOT MKL_CORE)
    set(MKL_LIBRARIES MKL_LIBRARIES-NOTFOUND)
    return()
endif()

get_filename_component(MKL_LIBDIR ${MKL_CORE} PATH)

find_library(MKL_SEQUENTIAL mkl_sequential ${MKL_LIBDIR})

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)

    if (APPLE)
        set(MKL_LIBRARIES MKL_LIBRARIES-NOTFOUND)
        return()
    endif()

	find_library(MKL_LP64 mkl_gf_lp64      ${MKL_LIBDIR})
	find_library(MKL_THREAD mkl_gnu_thread ${MKL_LIBDIR})
	list(APPEND MKL_THREAD -lgomp;-pthread)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI OR CMAKE_C_COMPILER_ID STREQUAL PGI)

	find_library(MKL_LP64 mkl_intel_lp64   ${MKL_LIBDIR})
	find_library(MKL_THREAD mkl_pgi_thread ${MKL_LIBDIR})

endif()


set(MKL_LIBRARIES ${MKL_LP64} ${MKL_THREAD} ${MKL_CORE})
set(LAPACK_LIBRARIES ${MKL_LP64} ${MKL_SEQUENTIAL} ${MKL_CORE})

message(STATUS "Found Intel MKL LP64 in ${MKL_LIBDIR}")
