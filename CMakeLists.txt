# CMakeLists.txt -- Build system for the Layer Lab
#
# Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.

cmake_minimum_required(VERSION 2.8)

project(layerlab)

option(USE_AVX "Enable AVX optimizations" TRUE)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE MinSizeRel CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6)
find_package(PythonLibs 3 REQUIRED)
find_package(PythonInterp 3 QUIET)

# Build TBB
if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()

set(TBB_BUILD_SHARED          OFF)
set(TBB_BUILD_STATIC          ON)
set(TBB_BUILD_TBBMALLOC       OFF)
set(TBB_BUILD_TBBMALLOC_PROXY OFF)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ext/tbb" tbb)

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)
if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unsequenced")
  endif()
  if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
  endif()
endif()

# Search for FFTW3
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(FFTW)

if (PKG_FFTW_FOUND)
  add_definitions("-DHAVE_FFTW")
endif()

# Compile with AVX
if (USE_AVX)
  if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
  endif()
endif()

# Compile with compiler warnings turned on
if(MSVC)
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# Set platform-specific flags
if (WIN32)
  add_definitions("-D__WINDOWS__")
elseif(UNIX)
  if(APPLE)
    add_definitions("-D__OSX__")
  else()
    add_definitions("-D__LINUX__")
  endif()
endif()

include_directories(
  ${PYTHON_INCLUDE_DIR}
  ext/eigen
  ext/filesystem
  ext/pybind11/include
  ext/tbb/include
  ext/tinyformat
  include
)

add_library(layerlab SHARED
  # Layer lab files
  include/layer/color.h
  include/layer/math.h       src/math.cpp
  include/layer/hg.h         src/hg.cpp
  include/layer/fourier.h    src/fourier.cpp
  include/layer/microfacet.h src/microfacet.cpp
  include/layer/layer.h      src/layer.cpp
  include/layer/fresnel.h    src/fresnel.cpp
  include/layer/mmap.h       src/mmap.cpp
  include/layer/log.h        src/log.cpp
  include/layer/storage.h    src/storage.cpp
  include/layer/simd.h       src/simd.cpp

  # Python API
  src/python.cpp src/py_vector.cpp src/py_math.cpp src/py_fourier.cpp
  src/py_spline.cpp src/py_layer.cpp src/py_quad.cpp src/py_fresnel.cpp
  src/py_filesystem.cpp
)

add_custom_target(mkdoc COMMAND
  python3 ${CMAKE_CURRENT_SOURCE_DIR}/ext/pybind11/tools/mkdoc.py 
  -I${PYTHON_INCLUDE_DIR} -Iext/eigen -Iext/pybind11/include -Iext/filesystem
  -Iext/tinyformat -Iinclude
  ${CMAKE_CURRENT_SOURCE_DIR}/include/layer/*.h > ${CMAKE_CURRENT_SOURCE_DIR}/src/py_doc.h)

set_target_properties(layerlab PROPERTIES PREFIX "")
link_directories(${CMAKE_CURRENT_BUILD_DIR}/tbb)
target_link_libraries(layerlab tbb_static)

if (PKG_FFTW_FOUND)
  target_link_libraries(layerlab ${FFTW_THREADS_LIB})
endif()

if (WIN32)
  if (MSVC)
    # Enforce size-based optimization and link time code generation on MSVC (~30% smaller binaries in experiments)
    set_target_properties(layerlab PROPERTIES COMPILE_FLAGS "/Os /GL")
    set_target_properties(layerlab PROPERTIES LINK_FLAGS "/LTCG")
  endif()

  # .PYD file extension on Windows
  set_target_properties(layerlab PROPERTIES SUFFIX ".pyd")

  # Link against the Python shared library
  target_link_libraries(layerlab ${PYTHON_LIBRARY})
elseif (UNIX)
  # It's quite common to have multiple copies of the same Python version
  # installed on one's system. E.g.: one copy from the OS and another copy
  # that's statically linked into an application like Blender or Maya.
  # If we link our plugin library against the OS Python here and import it
  # into Blender or Maya later on, this will cause segfaults when multiple
  # conflicting Python instances are active at the same time.

  # Windows does not seem to be affected by this issue. The solution for Linux
  # and Mac OS is simple: we just don't link against the Python library. The
  # resulting shared library will have missing symbols, but that's perfectly
  # fine -- they will be resolved at import time.

  # .SO file extension on Linux/Mac OS
  set_target_properties(layerlab PROPERTIES SUFFIX ".so")

  # Strip unnecessary sections of the binary on Linux/Mac OS
  if(APPLE)
    set_target_properties(layerlab PROPERTIES MACOSX_RPATH ".")
    set_target_properties(layerlab PROPERTIES LINK_FLAGS "-undefined dynamic_lookup -dead_strip")
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET layerlab POST_BUILD COMMAND strip -u -r layerlab.so)
    endif()
  else()
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET layerlab POST_BUILD COMMAND strip layerlab.so)
    endif()
  endif()
endif()
