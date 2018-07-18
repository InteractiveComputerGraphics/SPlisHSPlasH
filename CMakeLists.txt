cmake_minimum_required(VERSION 2.8)

project(SPlishSPlasH)

set(PROJECT_PATH ${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(${PROJECT_PATH}/CMake/Common.cmake)

if (WIN32)
  subdirs(
  extern/freeglut 
  extern/AntTweakBar 
  extern/glew 
  extern/zlib 
  extern/partio 
  extern/MD5
  Demos SPlisHSPlasH Utilities Tools)
else()
  subdirs(
	extern/AntTweakBar 
	extern/glew 
	extern/zlib 
	extern/partio 
	extern/md5
	Demos SPlisHSPlasH Utilities Tools)
endif()

add_definitions(-DSPH_DATA_PATH="../data")

OPTION(USE_DOUBLE_PRECISION "Use double precision"	ON)
if (USE_DOUBLE_PRECISION)
	add_definitions( -DUSE_DOUBLE)	
endif (USE_DOUBLE_PRECISION)

## PositionBasedDynamics
include(ExternalProject)
set(ExternalInstallDir "${CMAKE_SOURCE_DIR}/extern/install")
ExternalProject_Add(
   Ext_PBD
   PREFIX "${CMAKE_SOURCE_DIR}/extern/PositionBasedDynamics"
   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/PositionBasedDynamics.git
   GIT_TAG "b4d4c318947448bb13051a2b77713ecb58ad73b0"
   INSTALL_DIR ${ExternalInstallDir}/PositionBasedDynamics
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/PositionBasedDynamics -DPBD_NO_DEMOS:BOOL=1	-DPBD_EXTERNALINSTALLDIR:PATH=${ExternalInstallDir} -DUSE_DOUBLE_PRECISION:BOOL=${USE_DOUBLE_PRECISION}
) 

## CompactNSearch
ExternalProject_Add(
   Ext_CompactNSearch
   PREFIX "${CMAKE_SOURCE_DIR}/extern/CompactNSearch"
   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/CompactNSearch.git
   GIT_TAG "4ad505f5457e34ec479994740b75767241cff468"
   INSTALL_DIR ${ExternalInstallDir}/CompactNSearch
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/CompactNSearch -DUSE_DOUBLE_PRECISION:BOOL=${USE_DOUBLE_PRECISION}
) 

## GenericParameters
ExternalProject_Add(
   Ext_GenericParameters
   PREFIX "${CMAKE_SOURCE_DIR}/extern/GenericParameters"
   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/GenericParameters.git
   GIT_TAG "1ec904bf8555e78ae0d8ba2f9f9a395371c5d4eb"
   INSTALL_DIR ${ExternalInstallDir}/GenericParameters
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/GenericParameters -DGENERICPARAMETERS_NO_TESTS:BOOL=1		
) 

