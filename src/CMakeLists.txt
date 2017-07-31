cmake_minimum_required(VERSION 2.8)
 
project(ImageToImageFilter)

#---------------------------------------------
# Add FindTBB directory to CMake's module path
#---------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")
MESSAGE("CMAKE_SOURCE_DIR is ${CMAKE_SOURCE_DIR}")
MESSAGE("CMAKE_MODULE_PATH is ${CMAKE_MODULE_PATH}")

#---------------------------------------------
# Find required packages
#---------------------------------------------
FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(ITK REQUIRED)

IF(NOT TBB_FOUND)
        MESSAGE ( ERROR "TBB not found. Example project will not be built.")
ENDIF(NOT TBB_FOUND)

IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
  INCLUDE_DIRECTORIES(${USE_ITK_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

#---------------------------------------------
# Make it compatible with /trunk/ CRL code
#---------------------------------------------
SET(CRKIT_VERSION_MAJOR "1")
SET(CRKIT_VERSION_MINOR "5")
SET(CRKIT_VERSION_PATCH "0")
SET(CRKIT_VERSION_STRING "${CRKIT_VERSION_MAJOR}.${CRKIT_VERSION_MINOR} (${SVN_VERSION})")
CONFIGURE_FILE( configuration.h.in ${PROJECT_BINARY_DIR}/configuration.h)


INCLUDE_DIRECTORIES( libs/
                     libs/teem
                     libs/include 
                     libs/itkNLOPTOptimizers 
                     libs/crlDWICommon
                     libs/tclap
                     libs/eigen
                     ${CMAKE_SOURCE_DIR}/itk
                     ${CMAKE_SOURCE_DIR}/tbb
                     ${CMAKE_CURRENT_BINARY_DIR}
                     ${CMAKE_CURRENT_BINARY_DIR}/libs/itkNLOPTOptimizers )

#--------------------------------------------------------------------
# USE OPTIMIZED/NOT OPTIMIZED VERSION OF NLOPT
#--------------------------------------------------------------------
SET(NLOPTVERSION "original" CACHE STRING "Version of NLOPT")
SET(NLOPTVERSIONVALUES "original;optimized1;optimized2" CACHE STRING "List of possible values for NLOPT" )
SET_PROPERTY(CACHE NLOPTVERSION PROPERTY STRINGS ${NLOPTVERSIONVALUES}  )
MESSAGE(STATUS "NLOPTVERSION='${NLOPTVERSION}'")


SET(NLOPT_DIR nlopt-2.4.1-${NLOPTVERSION} )

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/libs/itkNLOPTOptimizers/${NLOPT_DIR} 
                     ${CMAKE_CURRENT_BINARY_DIR}/libs/itkNLOPTOptimizers/${NLOPT_DIR}  )


#---------------------------------------------
# The libs
#---------------------------------------------
#ADD_SUBDIRECTORY( libs/include  )
#ADD_SUBDIRECTORY( libs/itkNLOPTOptimizers  )
#ADD_SUBDIRECTORY( libs/crlDWICommon )
#ADD_SUBDIRECTORY( libs/teem  )

ADD_SUBDIRECTORY( "${CMAKE_SOURCE_DIR}/libs/include" )


SET(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/../../../install/examples)
#---------------------------------------------
# Set CMAKE RPATH
#---------------------------------------------
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 

SET(CMAKE_INSTALL_RPATH ${TBB_LIBRARY_DIRS}:${TBB_DEBUG_LIBRARY_DIRS})

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#---------------------------------------
# Define the source/header files
#---------------------------------------
SET(HDRS
    itk/itkRicianNoiseCorrectionFilter_itk.h
    itk/itkTensorReconstructionFilter_itk.h
    itk/itkPSDTensorEstimationFilter_itk.h

    itknd/itkRicianNoiseCorrectionFilter_itknd.h
    itknd/itkTensorReconstructionFilter_itknd.h
    itknd/itkPSDTensorEstimationFilter_itknd.h


    tbb/itkRicianNoiseCorrectionFilter_tbb.h
    tbb/itkTensorReconstructionFilter_tbb.h
    tbb/itkPSDTensorEstimationFilter_tbb.h

    tbbnd/itkRicianNoiseCorrectionFilter_tbbnd.h
    tbbnd/itkTensorReconstructionFilter_tbbnd.h
    tbbnd/itkPSDTensorEstimationFilter_tbbnd.h

    PositiveTensorFunctor.h

    itknd/itkCRLImageToImageFilter.h
    tbb/itkTBBImageToImageFilter.h
    tbbnd/itkTBBNDImageToImageFilter.h
)

ADD_EXECUTABLE(
    TestFilters
    itkTestFilters.cxx

    itk/itkRicianNoiseCorrectionFilter_itk.h
    itk/itkTensorReconstructionFilter_itk.h
    itk/itkPSDTensorEstimationFilter_itk.h
    itk/itkRicianNoiseCorrectionFilter_itk.txx
    itk/itkTensorReconstructionFilter_itk.txx
    itk/itkPSDTensorEstimationFilter_itk.txx

    itknd/itkRicianNoiseCorrectionFilter_itknd.h
    itknd/itkTensorReconstructionFilter_itknd.h
    itknd/itkPSDTensorEstimationFilter_itknd.h
    itknd/itkRicianNoiseCorrectionFilter_itknd.txx
    itknd/itkTensorReconstructionFilter_itknd.txx
    itknd/itkPSDTensorEstimationFilter_itknd.txx

    tbb/itkRicianNoiseCorrectionFilter_tbb.h
    tbb/itkTensorReconstructionFilter_tbb.h
    tbb/itkPSDTensorEstimationFilter_tbb.h
    tbb/itkRicianNoiseCorrectionFilter_tbb.txx
    tbb/itkTensorReconstructionFilter_tbb.txx
    tbb/itkPSDTensorEstimationFilter_tbb.txx

    tbbnd/itkRicianNoiseCorrectionFilter_tbbnd.h
    tbbnd/itkTensorReconstructionFilter_tbbnd.h
    tbbnd/itkPSDTensorEstimationFilter_tbbnd.h
    tbbnd/itkRicianNoiseCorrectionFilter_tbbnd.txx
    tbbnd/itkTensorReconstructionFilter_tbbnd.txx
    tbbnd/itkPSDTensorEstimationFilter_tbbnd.txx

    itknd/itkCRLImageToImageFilter.h
    itknd/itkCRLImageToImageFilter.txx
    tbb/itkTBBImageToImageFilter.h
    tbb/itkTBBImageToImageFilter.txx
    tbbnd/itkTBBNDImageToImageFilter.h
    tbbnd/itkTBBNDImageToImageFilter.txx
)

TARGET_LINK_LIBRARIES( TestFilters
    ${ITK_LIBRARIES}  
    ${TBB_LIBRARIES})

INSTALL(TARGETS  TestFilters
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)


