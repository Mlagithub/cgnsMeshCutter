
add_subdirectory(cgio)

# metis
set(METIS_ROOT "/home/one/.local/metis-5.1.0")
find_library(metisLib metis "${METIS_ROOT}/lib" NO_DEFAULT_PATH)
message(STATUS "Found metis: " ${metisLib})
target_link_libraries(${TgtLib} PUBLIC ${metisLib})
target_include_directories(${TgtLib} PUBLIC "${METIS_ROOT}/include")

# parmetis
set(PARMETIS_ROOT "/home/one/.local/parmetis-4.0.3")
find_library(parmetisLib parmetis "${PARMETIS_ROOT}/lib" NO_DEFAULT_PATH)
message(STATUS "Found parmetis: " ${parmetisLib})
target_link_libraries(${TgtLib} PUBLIC ${parmetisLib})
target_include_directories(${TgtLib} PUBLIC "${PARMETIS_ROOT}/include")

target_include_directories(${TgtLib} 
PUBLIC 
    ${CMAKE_CURRENT_LIST_DIR}
)

target_sources(${TgtLib}
PRIVATE
    meshCutter.cpp
    timeStat.cpp
    parmetisCutter.cpp
    mpiAdapter.cpp
    cartesianCutter.cpp
    stringUtil.cpp
)