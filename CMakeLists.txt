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


## PositionBasedDynamics
include(ExternalProject)
set(ExternalInstallDir "${CMAKE_SOURCE_DIR}/extern/install")
ExternalProject_Add(
   ExternalProject_PositionBasedDynamics
   PREFIX "${CMAKE_SOURCE_DIR}/extern/PositionBasedDynamics"
   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/PositionBasedDynamics.git
   GIT_TAG "6315c0ef206b028afc092ed191bb2c717ccd15bf"
   INSTALL_DIR ${ExternalInstallDir}/PositionBasedDynamics
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/PositionBasedDynamics -DPBD_NO_DEMOS:BOOL=1	-DPBD_EXTERNALINSTALLDIR:PATH=${ExternalInstallDir}
) 

## CompactNSearch
ExternalProject_Add(
   ExternalProject_CompactNSearch
   PREFIX "${CMAKE_SOURCE_DIR}/extern/CompactNSearch"
   GIT_REPOSITORY https://www.graphics.rwth-aachen.de:9000/jbender/CompactNSearch.git
   GIT_TAG "efc0f6a561ea18e903f1e57e13f6d5dcbc587dc5"
   INSTALL_DIR ${ExternalInstallDir}/CompactNSearch
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/CompactNSearch
) 

## GenericParameters
ExternalProject_Add(
   ExternalProject_GenericParameters
   PREFIX "${CMAKE_SOURCE_DIR}/extern/GenericParameters"
   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/GenericParameters.git
   GIT_TAG "1ec904bf8555e78ae0d8ba2f9f9a395371c5d4eb"
   INSTALL_DIR ${ExternalInstallDir}/GenericParameters
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/GenericParameters -DGENERICPARAMETERS_NO_TESTS:BOOL=1		
) 

