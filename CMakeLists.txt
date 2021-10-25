cmake_minimum_required(VERSION 3.0.0)
project(CGNSMeshCutter VERSION 0.1.0)

cmake_policy(SET CMP0076 NEW)
cmake_policy(SET CMP0079 NEW)

# 
set(TgtLib "cgnsCutter")
add_library(${TgtLib} SHARED "")

# Cutter
set(TgtExe "cutter.exe")
add_executable(${TgtExe} "")
target_link_libraries(${TgtExe} PUBLIC ${TgtLib})

add_subdirectory(src)
add_subdirectory(test)
