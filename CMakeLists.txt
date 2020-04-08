################################################################################
# general CMake and project setup
################################################################################
cmake_minimum_required(VERSION 3.1)

project(SPlishSPlasH)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(Common)
if (NOT WIN32)
	message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
endif()

set(ExternalInstallDir "${CMAKE_BINARY_DIR}/extern/install" CACHE INTERNAL "")
set(EXT_CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE INTERNAL "")
if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	set(EXT_CMAKE_BUILD_TYPE "Release" CACHE INTERNAL "")
endif()

set(TOPLEVEL_INCLUDE_DIR ${PROJECT_SOURCE_DIR})

################################################################################
# foreign external libraries
################################################################################
add_subdirectory(extern/zlib)
add_subdirectory(extern/partio)
add_subdirectory(extern/md5)
add_subdirectory(extern/tinyexpr)
if (NOT SPH_LIBS_ONLY)
	add_subdirectory(extern/AntTweakBar)
	add_subdirectory(extern/glew)
	if(WIN32)
		add_subdirectory(extern/freeglut)
	endif()
	if (USE_PYTHON_BINDINGS)
		add_subdirectory(extern/pybind)
	endif ()
endif()

## Eigen3 is used by most of the libraries that follow
find_package(Eigen3 REQUIRED)

################################################################################
# own external libraries
################################################################################
include(SetUpExternalProjects)

################################################################################
# internal libraries
################################################################################
add_subdirectory(SPlisHSPlasH)
add_subdirectory(Utilities)

################################################################################
# executables
################################################################################
if (NOT SPH_LIBS_ONLY)
	include(DataCopyTargets)
	add_subdirectory(Tools)
	add_subdirectory(Simulator)
	add_subdirectory(Tests)
	if (USE_PYTHON_BINDINGS)
		add_subdirectory(pySPlisHSPlasH)
	endif ()
endif()
