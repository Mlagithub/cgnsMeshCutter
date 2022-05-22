target_include_directories(${TgtLib} 
PUBLIC 
    ${CMAKE_CURRENT_LIST_DIR}
)

# CGNS
set(CGNS_ROOT "/home/one/.local/cgns-develop/")
find_library(cgnsLib cgns "${CGNS_ROOT}/lib" NO_DEFAULT_PATH)
message(STATUS "Found cgns: " ${cgnsLib})
target_link_libraries(${TgtLib} PUBLIC ${cgnsLib})
target_include_directories(${TgtLib} PUBLIC "${CGNS_ROOT}/include")

# fmt
# find_package(fmt CONFIG REQUIRED)
# if(fmt_FOUND)
#     message(STATUS "Found fmt: " ${fmt_DIR}/../../libfmt.so)
#     target_link_libraries(${TgtLib} PUBLIC ${fmt_DIR}/../../libfmt.so)
#     target_include_directories(${TgtLib} PUBLIC ${fmt_DIR}/../../../include) 
# endif(fmt_FOUND)

# MPI
find_package(MPI REQUIRED)
if(MPI_FOUND)
    message(STATUS "Found mpi: " ${MPI_LIBRARY})
    target_link_libraries(${TgtLib} PUBLIC ${MPI_LIBRARY})
    target_include_directories(${TgtLib} PUBLIC ${MPI_INCLUDE_PATH})
endif(MPI_FOUND)


target_sources(${TgtLib}
PRIVATE
    CGError.cpp
    CGFile.cpp
    CGIONode.cpp 
    CGNSBase.cpp 
    CGNSNode.cpp 
    Section.cpp 
    Zone.cpp 
)