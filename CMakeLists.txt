CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)
PROJECT(
   pgvl
   VERSION 0.1.0
   LANGUAGES CXX C
)

#==============Set Directories=================
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" )

#==============Custom Cmake Modules=============
#SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules" )

#================Find Packages==================
FIND_PACKAGE( GTest REQUIRED )
INCLUDE_DIRECTORIES( ${GTEST_INCLUDE_DIR} )
ENABLE_TESTING()

#=================Options=======================

OPTION(
   DO_PROFILING
   "Build with profiling flags set"
   OFF
)

#===========Set release/debug settings==========

SET( CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall -Wno-sign-compare" )
SET( LINK_FLAGS_DEBUG "-rdynamic" )
SET( CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wno-deprecated -Wno-sign-compare" )
SET( LINK_FLAGS_RELEASE "-rdynamic" )
IF( ${DO_PROFILING} )
   SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg" )
   SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pg" )
ENDIF()

# C++11
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

# OpenMP
SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )

IF( CMAKE_BUILD_TYPE STREQUAL "Release" )
  # This definition disables the assert() macro.
  ADD_DEFINITIONS( -DNDEBUG )
ENDIF()

#==============Create config.h==================
CONFIGURE_FILE( config.in config.h )
# config.h will be put in the binary dir, so have to add it as an include dir.
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

#============Add Include directories============
INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/include" )

#=============Process Subdirectories============
ADD_SUBDIRECTORY( "doc/" )
ADD_SUBDIRECTORY( "src/" )
ADD_SUBDIRECTORY( "tests/" )