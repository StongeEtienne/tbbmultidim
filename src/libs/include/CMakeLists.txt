# DEFINE PROJECT VARIABLES
# - provides CRLCOMMON_BINARY_DIR CRLCOMMON_SOURCE_DIR
PROJECT(CRLCOMMON)
 
FIND_PACKAGE(ITK)
IF (NOT ITK_FOUND)
  MESSAGE(FATAL_ERROR "Please set ITK_DIR.")
ENDIF(NOT ITK_FOUND)
INCLUDE(${ITK_USE_FILE})
INCLUDE_DIRECTORIES(${USE_ITK_FILE})


INCLUDE_DIRECTORIES( ../itkNLOPTOptimizers )

# Some classes depend on the compile time generated file: configuration.h
# which is available in ${CRKIT_BINARY_DIR}.
INCLUDE_DIRECTORIES(
${CRKIT_SOURCE_DIR}/common/external
${CRKIT_SOURCE_DIR}/common/include
${CRKIT_BINARY_DIR}
)

#---------------------------------------
# Define the source/header files
#---------------------------------------
SET(CRLCOMMON_HDRS
    crlFileName.h
    crlFileUtil.h
    crlLogClass.h
    crlGeneralizedInstantiator.h
    crlGeneralizedInstantiatorMacros.h
    crlNHDRFileUtils.h
    crlTaskProgress.h
    crlTaskProgress2.h
    crlVectorImageUtils.h
    crlTransformUtils.h
    crlDocumentedToolOutput.h
    crlRobustMeanCostFunction.h
    crlRobustMeanImageFilter.h
    crlFastOps.h
    crlMFMUtils.h
    crlOptimizerUtils.h   
    crlDilateLabelImageFilter.h
    crlDWIStudy.h
    crlDWIStudy.txx 
    crlDWICommon.h
    crlDWIGradients.h
    crlDWIGradientDirectionsDef.h
    crlDWIUtils.h
   )

SET(CRLCOMMON_SRCS
    crlDWICommon.cxx
    crlRobustMeanCostFunction.cxx
    crlDWIGradients.cxx
    crlFileName.cxx
    crlFileUtil.cxx
    crlLogClass.cxx
    crlNHDRFileUtils.cxx
    crlTaskProgress.cxx
    crlTaskProgress2.cxx
    crlFastOps.cxx
    crlOptimizerUtils.cxx    
    crlDilateLabelImageFilter.txx
    crlRobustMeanImageFilter.txx
   )

#---------------------------------------
# ONLY if ITK >= 4
#---------------------------------------
IF( "${ITK_VERSION_MAJOR}" GREATER 3 )

  SET(CRLCOMMON_SRCS
      ${CRLCOMMON_SRCS}
      crlMFMIO.cxx )

  SET(CRLCOMMON_HDRS
      ${CRLCOMMON_HDRS}
      crlMFMIO.h )

ENDIF()

#---------------------------------------
# Create the Doxygen  documentation
# (the path of the files must be absolute)
#---------------------------------------
OPTION(CRKIT_COMPILE_DOC "Compile the CRKit classes documentation with Doxygen" OFF )

IF(CRKIT_COMPILE_DOC)
  SET(TMP ${CRLCOMMON_HDRS} ${CRLCOMMON_SRCS} )
  SET(DOXYGEN_INPUT "" )
  FOREACH(I ${TMP})
     SET(DOXYGEN_INPUT "${DOXYGEN_INPUT} ${CMAKE_CURRENT_SOURCE_DIR}/${I}" )
  ENDFOREACH(I)

  SET( DOXYGEN_EXCLUDE  )

  SET( DOXYGEN_TARGET_NAME "doc_crlCommon" )
  SET( DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/docs/crlCommon")
  SET( DOXYGEN_INSTALL_DIR "${PROJECT_BINARY_DIR}/docs/crlCommon" )
  INCLUDE(../../DoxygenDoc.cmake)
ENDIF(CRKIT_COMPILE_DOC)

#---------------------------------------
# Add the library
#---------------------------------------
ADD_LIBRARY (crlCommon STATIC
  ${CRLCOMMON_SRCS}
  ${CRLCOMMON_HDRS}
)
LINK_DIRECTORIES(${CRLCOMMON_BINARY_DIR})
INSTALL(TARGETS crlCommon DESTINATION bin)

