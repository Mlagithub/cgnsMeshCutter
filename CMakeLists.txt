cmake_minimum_required(VERSION 3.0.0)
project(CGNSMeshCutter VERSION 0.1.0)

cmake_policy(SET CMP0076 NEW)
cmake_policy(SET CMP0079 NEW)

# compiler 
set(CMAKE_CXX_STANDARD 11)
set(buildType "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${buildType}'.")
    set(CMAKE_BUILD_TYPE "${buildType}" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# 
set(TgtLib "cgnsCutter")
add_library(${TgtLib} SHARED "")

# Cutter
set(TgtExe "cutter.exe")
add_executable(${TgtExe} "")
target_link_libraries(${TgtExe} PUBLIC ${TgtLib})

add_subdirectory(src)
add_subdirectory(test)
