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


# MPI
find_package(MPI)
if(MPI_FOUND)
    message(STATUS "Found MPI_LIB: " ${MPI_CXX_LIBRARIES})
    target_include_directories(${TgtLib} PUBLIC ${MPI_C_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_DIRS})
    target_link_libraries(${TgtLib} PUBLIC ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
else(MPI_FOUND)
    set (MPI_ROOT "$ENV{MPI_ROOT}")
    find_library(mpilib NAMES mpi PATHS "${MPI_ROOT}/lib" `which mpicc`../../lib)
    target_include_directories(${TgtLib} PUBLIC "${MPI_ROOT}/include")
    target_link_libraries(${TgtLib} PUBLIC ${mpilib})
endif(MPI_FOUND)

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