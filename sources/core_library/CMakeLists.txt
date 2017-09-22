cmake_minimum_required(VERSION 2.8)

project(LibTri-IA2D-core)
set(LibTri-IA2D_VERSION_MAJOR 1)
set(LibTri-IA2D_VERSION_MINOR 0)

set(LANGUAGE C++)
#set(CMAKE_VERBOSE_MAKEFILE on)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

set(CMAKE_BUILD_TYPE Release)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11") ## Optimize
#    set(CMAKE_EXE_LINKER_FLAGS "-s -lrt") ## Strip binary
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-O3 -march=native -std=c++11") ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s") ## Strip binary
endif()

include_directories(sources)
add_subdirectory(sources)

target_include_directories(_IA2D PUBLIC ${PROJECT_BINARY_DIR}/sources/)
