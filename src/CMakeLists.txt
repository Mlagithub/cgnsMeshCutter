# CGNS
set (CGNS_ROOT "$ENV{CGNS_ROOT}")
find_library(cgnslib cgns "${CGNS_ROOT}/lib" REQUIRED NO_DEFAULT_PATH)
if(cgnslib)
    message(STATUS "Found cgns: " ${cgnslib})
    target_link_libraries(${TgtLib} PUBLIC ${cgnslib})
    target_include_directories(${TgtLib} PUBLIC "${CGNS_ROOT}/include")
endif()

# METIS
set (METIS_ROOT "$ENV{METIS_ROOT}")
find_library(metislib NAMES metis PATHS "${METIS_ROOT}/lib" NO_DEFAULT_PATH)
if(metislib)
    message(STATUS "Found metis: " ${metislib})
    target_compile_definitions(${TgtLib} PUBLIC "SUPPORT_METIS")
    target_link_libraries(${TgtLib} PUBLIC ${metislib})
    target_include_directories(${TgtLib} PUBLIC "${METIS_ROOT}/include")
endif()

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
)