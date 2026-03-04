#
# A small macro used for converting a list into a space
# separated string:
#
MACRO(TO_STRING _variable)
  SET(${_variable} "")
  FOREACH(_var  ${ARGN})
    SET(${_variable} "${${_variable}} ${_var}")
  ENDFOREACH()
  STRING(STRIP "${${_variable}}" ${_variable})
ENDMACRO()

#
# If 'variable' is empty it will be set to 'value'
#
MACRO(SET_IF_EMPTY _variable)
  IF("${${_variable}}" STREQUAL "")
    SET(${_variable} ${ARGN})
  ENDIF()
ENDMACRO()


MACRO(_FIND_LIBRARY _library_name)
  TO_STRING(_str ${ARGN})
  FIND_LIBRARY(${_library_name} ${ARGN})
  IF(${_library_name} MATCHES "-NOTFOUND")
    MESSAGE(STATUS "${_library_name} not found! The call was:")
    MESSAGE(STATUS "    FIND_LIBRARY(${_library_name} ${_str})")
  ELSE()
    MESSAGE(STATUS "Found ${${_library_name}}")
    target_link_libraries(${TgtLib} PUBLIC ${${_library_name}})
  ENDIF()
ENDMACRO()

MACRO(_FIND_PATH _path_name)
  TO_STRING(_str ${ARGN})
  FIND_PATH(${_path_name} ${ARGN})
  IF(${_path_name} MATCHES "-NOTFOUND")
    MESSAGE(STATUS "${_path_name} not found! The call was:")
    MESSAGE(STATUS "    FIND_PATH(${_path_name} ${_str})")
  ELSE()
    MESSAGE(STATUS "Found ${${_path_name}}")
    target_include_directories(${TgtLib} PUBLIC ${${_path_name}})
  ENDIF()
ENDMACRO()

# CGNS
SET_IF_EMPTY(CGNS_DIR "$ENV{CGNS_ROOT}")
_FIND_LIBRARY(CGNS_LIBRARY
    NAMES cgns
    HINTS ${CGNS_DIR}
    PATH_SUFFIXES lib64 lib
)
_FIND_PATH(CGNS_INCLUDE_DIR cgnslib.h
    HINTS ${CGNS_DIR}
    PATH_SUFFIXES include cgns include/cgns
)

# METIS
set (METIS_ROOT "$ENV{METIS_ROOT}")
find_library(metislib NAMES metis PATHS "${METIS_ROOT}/lib" NO_DEFAULT_PATH)
if(metislib)
    message(STATUS "Found metis: " ${metislib})
    target_compile_definitions(${TgtLib} PUBLIC "SUPPORT_METIS")
    target_link_libraries(${TgtLib} PUBLIC ${metislib})
    target_include_directories(${TgtLib} PUBLIC "${METIS_ROOT}/include")
elseif()
    message(ERROR "Can not find metis.")
endif()

# ParMETIS
set (PARMETIS_ROOT "$ENV{PARMETIS_ROOT}")
find_library(parmetislib NAMES parmetis PATHS "${PARMETIS_ROOT}/lib" REQUIRED NO_DEFAULT_PATH)
if(parmetislib)
    target_compile_definitions(${TgtLib} PUBLIC "SUPPORT_PARMETIS")
    message(STATUS "Found parmetis: " ${parmetislib})
    target_include_directories(${TgtLib} PUBLIC "${PARMETIS_ROOT}/include")
    target_link_libraries(${TgtLib} PUBLIC ${parmetislib})
elseif()
    message(ERROR "Can not find parmetis.")
endif()

# Link external libraries
target_link_libraries(${TgtLib} PUBLIC compi)

target_include_directories(${TgtLib} 
PUBLIC 
    ${CMAKE_CURRENT_LIST_DIR}
)

target_sources(${TgtLib}
PRIVATE
    meshCutter.cpp
    stringUtil.cpp
    cmdLine.cpp
    cgFile.cpp
    cartesianCutter.cpp
    metisCutter.cpp
    mpiAdapter.cpp
)
