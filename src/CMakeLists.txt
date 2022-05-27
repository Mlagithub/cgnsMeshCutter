set(CGNS_ROOT "/opt/cgns-develop/")
find_library(cgnsLib cgns "${CGNS_ROOT}/lib" NO_DEFAULT_PATH)
message(STATUS "Found cgns: " ${cgnsLib})
target_link_libraries(${TgtLib} PUBLIC ${cgnsLib})
target_include_directories(${TgtLib} PUBLIC "${CGNS_ROOT}/include")

target_include_directories(${TgtLib} 
PUBLIC 
    ${CMAKE_CURRENT_LIST_DIR}
)

target_sources(${TgtLib}
PRIVATE
    meshCutter.cpp
    cartesianCutter.cpp
    stringUtil.cpp
)