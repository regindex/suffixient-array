CMAKE_MINIMUM_REQUIRED(VERSION 3.13)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

# Set a default build type if none was specified
# ------------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# About this project
# ------------------------------------------------------------------------------
project (suffixient)
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "1")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

#option(COMPILE_BENCHMARKS "Set ON to compile the benchmarks" OFF) To be added what running experiments

# Set environment
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/install)

find_package(Git)
if(GIT_FOUND)
    message("git found: ${GIT_EXECUTABLE}")
else()
    message(WARNING "git not found. Cloning of submodules will not work.")
endif()

# Configure the compiler with the appropriate flags
# ------------------------------------------------------------------------------
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  # using Clang
  message("compiler: configuring Clang compiler")
  include(ConfigureCompilerClang)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  message("compiler: configuring gcc compiler")
  include(ConfigureCompilerGcc)
else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message(FATAL_ERROR "Only the compiler gcc and clang are supported")
endif()

# Configure thirdparty
# ------------------------------------------------------------------------------
set(CMAKE_INSTALL_INCLUDEDIR "include")
add_subdirectory(thirdparty)

add_subdirectory(include)
add_subdirectory(suff-set-src)
add_subdirectory(sA-index-src)
add_subdirectory(pfp-src)

#add_subdirectory(experiments)

# Set script containing the full pipeline
# ------------------------------------------------------------------------------
configure_file(${PROJECT_SOURCE_DIR}/pipeline/suffixient-set.py ${PROJECT_BINARY_DIR}/suffixient-set.py)
configure_file(${PROJECT_SOURCE_DIR}/pipeline/suffixient-array-index.py ${PROJECT_BINARY_DIR}/suffixient-array-index.py)