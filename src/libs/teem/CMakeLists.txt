#
# Teem: Tools to process and visualize scientific data and images
# Copyright (C) 2013, 2012, 2011, 2010, 2009  University of Chicago
# Copyright (C) 2008, 2007, 2006, 2005  Gordon Kindlmann
# Copyright (C) 2004, 2003, 2002, 2001, 2000, 1999, 1998  University of Utah
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# (LGPL) as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# The terms of redistributing and/or modifying this software also
# include exceptions to the LGPL that facilitate static linking.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#

# This is the root Teem CMakeLists file
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

# Cmake shouldn't define WIN32 when compiling
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

# As of CMake 2.6 Policies were introduced in order to provide a mechanism for
# adding backwards compatibility one feature at a time.

# Run cmake --help-policy CMP<num> to see documentation.
IF(COMMAND cmake_policy)
  # Library paths (/path/to/libmy.so not translated to -L/path/to -lmy)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

PROJECT(Teem)

#-----------------------------------------------------------------------------
# Teem version number.  Must reflect the values in teem/src/air/air.h
SET(Teem_VERSION_MAJOR "1")
SET(Teem_VERSION_MINOR "11")
SET(Teem_VERSION_PATCH "1")

# Version string should not include patch level.  The major.minor is
# enough to distinguish available features of the toolkit.
SET(Teem_VERSION_STRING "${Teem_VERSION_MAJOR}.${Teem_VERSION_MINOR}.${Teem_VERSION_PATCH}")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMake")

# We need ansi c-flags, especially on HP
SET(CMAKE_C_FLAGS "${CMAKE_ANSI_CFLAGS} ${CMAKE_C_FLAGS}")
SET(CMAKE_REQUIRED_FLAGS ${CMAKE_ANSI_CFLAGS})

#-----------------------------------------------------------------------------
# Test for some required system information.
INCLUDE (CMakeBackwardCompatibilityC)

#-----------------------------------------------------------------------------
# Output directories.
#SET (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
#SET (EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
#output directory for installing all include files.")
#MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
#SET(Teem_LIBRARY_PATH "${LIBRARY_OUTPUT_PATH}")
#SET(Teem_EXECUTABLE_PATH "${EXECUTABLE_OUTPUT_PATH}")
#SET(C_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})


FILE(MAKE_DIRECTORY ${Teem_BINARY_DIR}/include)
FILE(MAKE_DIRECTORY ${Teem_BINARY_DIR}/include/teem)

#-----------------------------------------------------------------------------
# Find platform-specific differences in the handling of IEEE floating point
# special values.
INCLUDE(${Teem_SOURCE_DIR}/CMake/TestQnanhibit.cmake)
TEST_QNANHIBIT(QNANHIBIT_VALUE ${Teem_SOURCE_DIR}/CMake)
IF(QNANHIBIT_VALUE)
  SET(QNANHIBIT 1 CACHE INTERNAL "The 22nd bit of 32-bit floating-point quiet NaN.")
ELSE(QNANHIBIT_VALUE)
  SET(QNANHIBIT 0 CACHE INTERNAL "The 22nd bit of 32-bit floating-point quiet NaN.")
ENDIF(QNANHIBIT_VALUE)


#-----------------------------------------------------------------------------
# Teem build configuration options.
ADD_DEFINITIONS(-DTEEM_STATIC)


SET(Teem_ZLIB_LIB "")
SET(Teem_PNG_LIB "")
SET(Teem_ZLIB OFF )
SET(Teem_PNG OFF )
SET(Teem_BZIP2 OFF )
SET(Teem_PTHREAD OFF)
SET(Teem_LEVMAR OFF )
SET(Teem_FFTW3 OFF )

#-----------------------------------------------------------------------------

IF(${CMAKE_GENERATOR} MATCHES "Visual Studio")
  # For Visual Studio we don't care about warnings about deprecated sprintf and
  # the like.
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio")

#Teem Defines
#The QNANHIBIT variable is configured by the root level CMakeLists.txt
IF(QNANHIBIT)
  ADD_DEFINITIONS(-DTEEM_QNANHIBIT=1)
ELSE(QNANHIBIT)
  ADD_DEFINITIONS(-DTEEM_QNANHIBIT=0)
ENDIF(QNANHIBIT)


ADD_DEFINITIONS(-DTEEM_DIO=0)


# (TEEM_LIB_LIST)
SET(Teem_LIBRARIES air hest biff nrrd ell unrrdu moss gage dye limn echo hoover seek ten pull mite meet)

SET(Teem_HEADER_FILES
  air/air.h
  hest/hest.h
  biff/biff.h
  nrrd/nrrd.h   nrrd/nrrdDefines.h   nrrd/nrrdMacros.h  nrrd/nrrdEnums.h
  ell/ell.h     ell/ellMacros.h
  unrrdu/unrrdu.h
  moss/moss.h
  gage/gage.h
  dye/dye.h
  limn/limn.h
  echo/echo.h
  hoover/hoover.h
  seek/seek.h
  ten/ten.h     ten/tenMacros.h
  pull/pull.h
  mite/mite.h
  meet/meet.h
  )

SET(Teem_INSTALLED_HEADER_FILES)
FOREACH(header_file ${Teem_HEADER_FILES})
  # When debugging, uncomment this line
  #MESSAGE(STATUS "Copy header file: ${header_file}")
  GET_FILENAME_COMPONENT(file_no_path "${header_file}" NAME)
  CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/src/${header_file}"
    "${CMAKE_CURRENT_BINARY_DIR}/include/teem/${file_no_path}" COPYONLY IMMEDIATE)
  SET(Teem_INSTALLED_HEADER_FILES
    ${Teem_INSTALLED_HEADER_FILES}
    "${CMAKE_CURRENT_BINARY_DIR}/include/teem/${file_no_path}")
ENDFOREACH(header_file)

#---
# Include directory
INCLUDE_DIRECTORIES(
  "${Teem_BINARY_DIR}/include/"
  "${Teem_SOURCE_DIR}/include/"
  )


# All the source files
SET(Teem_SOURCES)

#-----------------------------------------------------------------------------
# Macro for adding the current set of source files to the master list
MACRO(ADD_Teem_LIBRARY dir)
  SET(srcs)
  FOREACH(src ${ARGN})
    SET(srcs ${srcs} "src/${dir}/${src}")
  ENDFOREACH(src ${ARGN})
  SET(Teem_SOURCES ${Teem_SOURCES} ${srcs})
  # This will group all the source files in the VS project by directory
  SOURCE_GROUP( ${dir} FILES ${srcs} )
ENDMACRO(ADD_Teem_LIBRARY name)


#-----------------------------------------------------------------------------
# Dispatch the build into the proper subdirectories.
FOREACH(dir ${Teem_LIBRARIES})
  INCLUDE(src/${dir}/sources.cmake)
ENDFOREACH(dir ${Teem_LIBRARIES})

#-----------------------------------------------------------------------------
# Build the Mega library
#
# This needs to occur after the parsing of the subdirectories.

# Create the library target
ADD_LIBRARY(teem ${Teem_SOURCES})

# Set up some library paths for installation.  Windows will ingore the
# RPATH stuff, and mac will ingore the INSTALL_RPATH.  Make sure for
# macs to set BUILD_WITH_INSTALL_RPATH OFF and set INSTALL_NAME_DIR.
# I belive INSTALL_NAME_DIR will be ignored on linux.
SET_TARGET_PROPERTIES(teem PROPERTIES
  BUILD_WITH_INSTALL_RPATH OFF
  INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib
  INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib
  SOVERSION 1
  VERSION ${Teem_VERSION_STRING}
  )


IF(UNIX)
  TARGET_LINK_LIBRARIES(teem -lm)
ENDIF(UNIX)


#-----------------------------------------------------------------------------
# Help outside projects build Teem projects.
INCLUDE(CMakeExportBuildSettings)
EXPORT_LIBRARY_DEPENDENCIES(${Teem_BINARY_DIR}/TeemLibraryDepends.cmake)
CMAKE_EXPORT_BUILD_SETTINGS(${Teem_BINARY_DIR}/TeemBuildSettings.cmake)

SET(CFLAGS "${CMAKE_C_FLAGS}")
SET(CC "${CMAKE_C_COMPILER}")

CONFIGURE_FILE("${Teem_SOURCE_DIR}/CMake/TeemConfig.cmake.in"
  "${Teem_BINARY_DIR}/TeemConfig.cmake" @ONLY IMMEDIATE)

