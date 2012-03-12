include (FindLibraryList)

if($ENV{BLA_VENDOR} MATCHES ".+")
  set(BLA_VENDOR $ENV{BLA_VENDOR} CACHE STRING "The BLAS Vendor" FORCE)
  message(${BLA_VENDOR} $ENV{BLA_VENDOR})
elseif(NOT BLA_VENDOR)
  set(BLA_VENDOR "All")
endif()

if (BLA_VENDOR MATCHES "ACML.*" OR BLA_VENDOR STREQUAL "All")
  if(NOT BLAS_LIBRARIES)
    # try to find acml in "standard" paths
    if ($ENV{ACML_ROOT} MATCHES ".+")
      set(_ACML_ROOT $ENV{ACML_ROOT})
    elseif(ACML_ROOT)
      set(_ACML_ROOT ${ACML_ROOT})
    else()
      file( GLOB _ACML_ROOT "/opt/acml*/ACML-EULA.txt" )

      if( _ACML_ROOT )
	list(GET _ACML_ROOT 0 _ACML_ROOT)
	get_filename_component( _ACML_ROOT ${_ACML_ROOT} PATH )
      endif()
    endif ()
    if( _ACML_ROOT )
      if( SIZEOF_INTEGER EQUAL 8 )
	set( _ACML_PATH_SUFFIX "_int64" )
      else()
	set( _ACML_PATH_SUFFIX "" )
      endif()
      if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
	set( _ACML_COMPILER "ifort64" )
      elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro" )
	set( _ACML_COMPILER "sun64" )
      elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
	set( _ACML_COMPILER "pgi64" )
      elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Open64" )
	set( _ACML_COMPILER "open64_64" )
      elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "NAG" )
	set( _ACML_COMPILER "nag64" )
      else() #if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
	set( _ACML_COMPILER "gfortran64" )
      endif()

      if( BLA_VENDOR MATCHES "ACML(_.*)" )
	string(TOLOWER ${CMAKE_MATCH_1} _ACML_VARIANT)
      endif()
      set(_lib "acml")
      if( BLA_VENDOR MATCHES "ACML.*MP" )
	set(_lib "acml_mp")
      endif()

      set(_ACML_LIBDIR "${_ACML_ROOT}/${_ACML_COMPILER}${_ACML_VARIANT}${_ACML_PATH_SUFFIX}/lib" )
      find_library_list(BLAS_LIBRARIES "${_lib}" ${_ACML_LIBDIR})
    endif()
  endif()
endif()

# Apple BLAS library?
if (BLA_VENDOR STREQUAL "Apple" OR BLA_VENDOR STREQUAL "All")
  if(NOT BLAS_LIBRARIES)
    find_library(BLAS_LIBRARIES "Accelerate")
  endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "Apple" OR BLA_VENDOR STREQUAL "All")

if (BLA_VENDOR MATCHES "mkl_.*" OR BLA_VENDOR STREQUAL "All")
  if( NOT BLAS_LIBRARIES)

    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      set(_MKL_LP64 "mkl_intel_lp64")
      set(_MKL_THREAD "mkl_intel_thread;mkl_core;iomp5")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(_MKL_LP64 "mkl_gf_lp64")
      set(_MKL_THREAD "mkl_gnu_thread;mkl_core;gomp")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI OR CMAKE_C_COMPILER_ID STREQUAL PGI)
      set(_MKL_LP64 "mkl_intel_lp64")
      set(_MKL_THREAD "mkl_pgi_thread;mkl_core")
    endif()

    if (BLA_VENDOR STREQUAL "mkl_sequential")
      set(_MKL_THREAD "mkl_sequential;mkl_core")
    endif()

    find_library_list(BLAS_LIBRARIES "${_MKL_LP64};${_MKL_THREAD}")
  endif()
endif()
