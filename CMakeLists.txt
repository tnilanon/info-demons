CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

PROJECT(InfoDemonWei)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
		"ITK NOT FOUND. PLEASE SET ITK_DIR.")
ENDIF(ITK_FOUND)


#if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-cpp")
#endif()



ADD_EXECUTABLE(InfoDemonMain MultiChannelDemon.cxx LabelDiffeomorphicDemonsRegistrationFilter.txx LabelDiffeomorphicDemonsRegistrationFilter.h LabelESMDemonsRegistrationFunction.txx LabelESMDemonsRegistrationFunction.h LabelMultiResolutionPDEDeformableRegistration.txx LabelMultiResolutionPDEDeformableRegistration.h HistogramField.txx HistogramField.h)

ADD_EXECUTABLE(BuildHistogramMain BuildHistogram.cxx LabelDiffeomorphicDemonsRegistrationFilter.txx LabelDiffeomorphicDemonsRegistrationFilter.h LabelESMDemonsRegistrationFunction.txx LabelESMDemonsRegistrationFunction.h LabelMultiResolutionPDEDeformableRegistration.txx LabelMultiResolutionPDEDeformableRegistration.h HistogramField.txx HistogramField.h SoftHistogram.h)

ADD_EXECUTABLE(NormalizeHistogramMain NormalizeHistogram.cxx LabelDiffeomorphicDemonsRegistrationFilter.txx LabelDiffeomorphicDemonsRegistrationFilter.h LabelESMDemonsRegistrationFunction.txx LabelESMDemonsRegistrationFunction.h LabelMultiResolutionPDEDeformableRegistration.txx LabelMultiResolutionPDEDeformableRegistration.h HistogramField.txx HistogramField.h SoftHistogram.h)

TARGET_LINK_LIBRARIES(InfoDemonMain ITKCommon ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(BuildHistogramMain ITKCommon ITKIO ITKStatistics)
TARGET_LINK_LIBRARIES(NormalizeHistogramMain ITKCommon ITKIO ITKStatistics)