set(SIMULATION_LINK_LIBRARIES Utilities AntTweakBar glew partio zlib SPlisHSPlasH MD5 tinyexpr)
set(SIMULATION_DEPENDENCIES Utilities AntTweakBar glew partio zlib SPlisHSPlasH MD5 tinyexpr)

if(WIN32)
  set(SIMULATION_LINK_LIBRARIES freeglut opengl32.lib glu32.lib ${SIMULATION_LINK_LIBRARIES})
  set(SIMULATION_DEPENDENCIES freeglut ${SIMULATION_DEPENDENCIES})
else()
  set(OpenGL_GL_PREFERENCE GLVND)
  find_package(GLUT REQUIRED)
  find_package(OpenGL REQUIRED)

  include_directories( ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})

  set(SIMULATION_LINK_LIBRARIES
	${SIMULATION_LINK_LIBRARIES}
	${GLUT_LIBRARIES}
	${OPENGL_LIBRARIES}
  )
endif()

find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

############################################################
# NeighborhoodSearch
############################################################
include_directories(${NEIGHBORHOOD_SEARCH_INCLUDE_DIR})
set(SIMULATION_DEPENDENCIES ${SIMULATION_DEPENDENCIES} Ext_NeighborhoodSearch)
set(SIMULATION_LINK_LIBRARIES ${SIMULATION_LINK_LIBRARIES} ${NEIGHBORHOOD_SEARCH_LIBRARIES})

############################################################
# DiscreGrid
############################################################
include_directories(${PROJECT_PATH}/extern/install/Discregrid/include)
set(SIMULATION_DEPENDENCIES ${SIMULATION_DEPENDENCIES} Ext_Discregrid)
set(SIMULATION_LINK_LIBRARIES ${SIMULATION_LINK_LIBRARIES}
	optimized Discregrid
	debug Discregrid_d)
link_directories(${PROJECT_PATH}/extern/install/Discregrid/lib)

############################################################
# GenericParameters
############################################################
include_directories(${PROJECT_PATH}/extern/install/GenericParameters/include)
set(SIMULATION_DEPENDENCIES ${SIMULATION_DEPENDENCIES} Ext_GenericParameters)

############################################################
# PBD
############################################################
include_directories(${PROJECT_PATH}/extern/install/PositionBasedDynamics/include)
add_definitions(-DUSE_PBD)
add_definitions(-DPBD_DATA_PATH="${PROJECT_PATH}/extern/install/PositionBasedDynamics/include/data")


set(PBDWRAPPER_SOURCE_FILES	PositionBasedDynamicsWrapper/PBDWrapper.cpp)
set(PBDWRAPPER_HEADER_FILES
	PositionBasedDynamicsWrapper/PBDWrapper.h
	PositionBasedDynamicsWrapper/PBDRigidBody.h)

set(SIMULATION_LINK_LIBRARIES ${SIMULATION_LINK_LIBRARIES}
	optimized Simulation debug Simulation_d
	optimized PositionBasedDynamics debug PositionBasedDynamics_d
	optimized Utils debug Utils_d)

set(SIMULATION_DEPENDENCIES Ext_PBD ${SIMULATION_DEPENDENCIES})
link_directories(${PROJECT_PATH}/extern/install/PositionBasedDynamics/lib)

add_executable(DynamicBoundarySimulator
	main.cpp

	${PROJECT_PATH}/Simulators/Common/SimulatorBase.cpp
	${PROJECT_PATH}/Simulators/Common/SimulatorBase.h
	${PROJECT_PATH}/Simulators/Common/TweakBarParameters.cpp
	${PROJECT_PATH}/Simulators/Common/TweakBarParameters.h

	${VIS_FILES}

	${PBDWRAPPER_SOURCE_FILES}
	${PBDWRAPPER_HEADER_FILES}
)

if(DL_OUTPUT)
	add_definitions(-DDL_OUTPUT)
endif()

add_definitions(-DTW_NO_LIB_PRAGMA -DTW_STATIC)

include_directories(${PROJECT_PATH}/extern/freeglut/include)
include_directories(${PROJECT_PATH}/extern/glew/include)

set_target_properties(DynamicBoundarySimulator PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(DynamicBoundarySimulator PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(DynamicBoundarySimulator PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
add_dependencies(DynamicBoundarySimulator ${SIMULATION_DEPENDENCIES})
target_link_libraries(DynamicBoundarySimulator ${SIMULATION_LINK_LIBRARIES})
VIS_SOURCE_GROUPS()

source_group("Header Files\\PBD" FILES ${PBDWRAPPER_HEADER_FILES})
source_group("Source Files\\PBD" FILES ${PBDWRAPPER_SOURCE_FILES})

set_target_properties(DynamicBoundarySimulator PROPERTIES FOLDER "Simulators")

add_custom_command(TARGET DynamicBoundarySimulator PRE_BUILD
					COMMAND ${CMAKE_COMMAND} -E copy_directory
					${PROJECT_PATH}/data/shaders $<TARGET_FILE_DIR:DynamicBoundarySimulator>/resources/shaders)

add_custom_command(TARGET DynamicBoundarySimulator PRE_BUILD
					COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${PROJECT_PATH}/extern/install/PositionBasedDynamics/include/data/shaders $<TARGET_FILE_DIR:DynamicBoundarySimulator>/resources/pbd_shaders)
