# CMakeLists.txt -- Build system for the Layer Lab
#
# Copyright (c) 2015 Wenzel Jakob <wenzel@inf.ethz.ch>
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.

cmake_minimum_required(VERSION 2.8)

project(layerlab)

option(USE_AVX "Enable AVX optimizations" TRUE)
set(LAYERLAB_PYTHON_VERSION "" CACHE STRING "Python version to use for compiling layer lab")

include(CheckCXXCompilerFlag)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE MinSizeRel CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)

# Try to autodetect Python (can be overridden manually if needed)
set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6 3.7)
if (NOT ${LAYERLAB_PYTHON_VERSION} STREQUAL "")
  find_package(PythonLibs ${LAYERLAB_PYTHON_VERSION} EXACT)
  if (NOT PythonLibs_FOUND)
    find_package(PythonLibs ${LAYERLAB_PYTHON_VERSION} REQUIRED)
  endif()
else()
  find_package(PythonLibs REQUIRED)
endif()

# Build TBB
if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif()

option(TBB_BUILD_SHARED          "" OFF)
option(TBB_BUILD_STATIC          "" ON)
option(TBB_BUILD_TBBMALLOC       "" OFF)
option(TBB_BUILD_TBBMALLOC_PROXY "" OFF)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ext/tbb" tbb)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  CHECK_CXX_COMPILER_FLAG("-std=c++14" HAS_CPP14_FLAG)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_CPP11_FLAG)

  if (HAS_CPP14_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  elseif (HAS_CPP11_FLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  else()
    message(FATAL_ERROR "Unsupported compiler -- at least C++11 support is needed!")
  endif()

  # Enable link time optimization and set the default symbol
  # visibility to hidden (very important to obtain small binaries)
  if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    # Default symbol visibility
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")

    # Check for Link Time Optimization support
    CHECK_CXX_COMPILER_FLAG("-flto" HAS_LTO_FLAG)
    if (HAS_LTO_FLAG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    endif()
  endif()
endif()

# Search for FFTW3
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(FFTW)

if (PKG_FFTW_FOUND)
  add_definitions("-DHAVE_FFTW")
  include_directories(${FFTW_INCLUDES})
endif()

# Compile with AVX
if (USE_AVX)
  if (MSVC)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
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
  add_definitions(-D__WINDOWS__ -D_UNICODE)
elseif(UNIX)
  if(APPLE)
    add_definitions(-D__OSX__)
  else()
    add_definitions(-D__LINUX__)
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
  include/layer/common.h
  include/layer/color.h
  include/layer/frame.h
  include/layer/spline.h
  include/layer/quad.h
  include/layer/simd.h
  include/layer/vector.h
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
    set_target_properties(layerlab PROPERTIES COMPILE_FLAGS "/Os /GL /MP /bigobj")
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
