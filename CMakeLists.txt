cmake_minimum_required(VERSION 3.1)

project(SPlishSPlasH)

set(PROJECT_PATH ${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(${PROJECT_PATH}/CMake/Common.cmake)

add_definitions(-DSPH_DATA_PATH="../data")

if (NOT WIN32)
	message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()

OPTION(USE_DOUBLE_PRECISION "Use double precision"	ON)
if (USE_DOUBLE_PRECISION)
	add_definitions( -DUSE_DOUBLE)
endif (USE_DOUBLE_PRECISION)

set(ExternalInstallDir "${CMAKE_SOURCE_DIR}/extern/install" CACHE INTERNAL "")
set(EXT_CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE INTERNAL "")
if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	set(EXT_CMAKE_BUILD_TYPE "Release" CACHE INTERNAL "")
endif()

include(${PROJECT_PATH}/CMake/NeighborhoodSearch.cmake)

if (WIN32)
  subdirs(
  extern/zlib
  extern/partio
  extern/MD5
  extern/tinyexpr
  SPlisHSPlasH Utilities)
else()
  subdirs(
	extern/zlib
	extern/partio
	extern/md5
	extern/tinyexpr
	SPlisHSPlasH Utilities)
endif()

if (NOT SPH_LIBS_ONLY)
	if (WIN32)
		add_subdirectory(extern/freeglut)
	endif()
	add_subdirectory(extern/AntTweakBar)
	add_subdirectory(extern/glew)
	add_subdirectory(Simulators)
	add_subdirectory(Tools)
	add_subdirectory(Tests)
endif()

## PositionBasedDynamics
include(ExternalProject)
ExternalProject_Add(
   Ext_PBD
   PREFIX "${CMAKE_SOURCE_DIR}/extern/PositionBasedDynamics"
   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/PositionBasedDynamics.git
   GIT_TAG "dd185b0a0479248cb0c1e98294532b55cb0d287d"
   INSTALL_DIR ${ExternalInstallDir}/PositionBasedDynamics
   CMAKE_ARGS -DCMAKE_BUILD_TYPE=${EXT_CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/PositionBasedDynamics -DPBD_NO_DEMOS:BOOL=1	-DPBD_EXTERNALINSTALLDIR:PATH=${ExternalInstallDir} -DUSE_DOUBLE_PRECISION:BOOL=${USE_DOUBLE_PRECISION}
)

## GenericParameters
ExternalProject_Add(
   Ext_GenericParameters
   PREFIX "${CMAKE_SOURCE_DIR}/extern/GenericParameters"
   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/GenericParameters.git
   GIT_TAG "b1ad669fac8d106515f6aa8514a03598d5766a36"
   INSTALL_DIR ${ExternalInstallDir}/GenericParameters
   CMAKE_ARGS -DCMAKE_BUILD_TYPE=${EXT_CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/GenericParameters -DGENERICPARAMETERS_NO_TESTS:BOOL=1		
) 

## Discregrid
ExternalProject_Add(
   Ext_Discregrid
   PREFIX "${CMAKE_SOURCE_DIR}/extern/Discregrid"
   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/Discregrid.git
   GIT_TAG "c0fb5aeac4c8a83e9f37c720315f13a834409b81"
   INSTALL_DIR ${ExternalInstallDir}/Discregrid
   CMAKE_ARGS -DCMAKE_BUILD_TYPE:STRING=${EXT_CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/Discregrid -DBUILD_CMD_EXECUTABLE:BOOL=0 -DEIGEN3_INCLUDE_DIR:PATH=${EIGEN3_INCLUDE_DIR}
) 
