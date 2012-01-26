if (APPLE)
	set(_libdir ENV DYLD_LIBRARY_PATH)
else()
	set(_libdir ENV LD_LIBRARY_PATH)
endif()

find_library(CHOLMOD_LIBRARIES NAMES cholmod PATHS ${_libdir})
if (NOT CHOLMOD_LIBRARIES)
	return()
endif()

get_filename_component(CHOLMOD_LIBDIR ${CHOLMOD_LIBRARIES} PATH)
get_filename_component(CHOLMOD_HOME ${CHOLMOD_LIBDIR} PATH)


foreach(lib amd camd colamd ccolamd)
	find_library(${lib}_lib NAMES "${lib}" PATHS ${CHOLMOD_HOME}/lib NO_DEFAULT_PATH)
	set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${${lib}_lib})
endforeach()

find_library(metis_lib NAMES metis PATHS ${_libdir})
if (metis_lib)
	set(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARIES} ${metis_lib})
endif()

find_path(CHOLMOD_INC cholmod.h PATH_SUFFIXES ufsparse suitesparse HINTS ${CHOLMOD_HOME}/include)

message(STATUS "Found CHOLMOD in ${CHOLMOD_HOME}")
message(STATUS "\theaders in ${CHOLMOD_INC}")
