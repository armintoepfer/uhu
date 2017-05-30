########################################
# CMake build script for the UHU library
########################################

cmake_policy(SET CMP0048 NEW)
project(UHU VERSION 0.0.1 LANGUAGES CXX C)
cmake_minimum_required(VERSION 3.2)

set(ROOT_PROJECT_NAME ${PROJECT_NAME} CACHE STRING "root project name")

# Build type
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile RelWithDebInfo ReleaseWithAssert" FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Build-time options
option(UHU_build_bin    "Build binaries." ON)
option(UHU_use_ccache   "Build UHU using ccache, if available." ON)

# Main project paths
set(UHU_RootDir       ${UHU_SOURCE_DIR})
set(UHU_IncludeDir    ${UHU_RootDir}/include)
set(UHU_SourceDir     ${UHU_RootDir}/src)
set(UHU_SwigDir       ${UHU_RootDir}/swig)
set(UHU_TestsDir      ${UHU_RootDir}/tests)
set(UHU_ThirdPartyDir ${UHU_RootDir}/third-party)

# Project configuration
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake ${CMAKE_MODULE_PATH})

# Fixed order, do not sort or shuffle
include(uhu-ccache)
include(uhu-releasewithassert)
include(uhu-dependencies)
include(uhu-compilerflags)
include(uhu-gitsha1)

# Build library
add_subdirectory(${UHU_SourceDir})