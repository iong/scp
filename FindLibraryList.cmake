macro(Find_Library_List LIBRARIES _list)
  if (WIN32)
    set(_libdir ENV LIB)
  elseif (APPLE)
    set(_libdir ENV DYLD_LIBRARY_PATH)
  else ()
    set(_libdir  ENV LD_LIBRARY_PATH)
  endif ()
  
  set(_libdir HINTS ${ARGN} ${_libdir})

  foreach(_lib ${_list})
    find_library(${_lib}_LIBRARY NAMES ${_lib} ${_libdir})
    set(${LIBRARIES} ${${LIBRARIES}} ${${_lib}_LIBRARY})
  endforeach()
endmacro()


