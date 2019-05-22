set(SIMULATION_LINK_LIBRARIES AntTweakBar glew SPlisHSPlasH Utilities partio zlib MD5)
set(SIMULATION_DEPENDENCIES AntTweakBar glew SPlisHSPlasH Utilities partio zlib MD5)

if(WIN32)
  set(SIMULATION_LINK_LIBRARIES freeglut opengl32.lib glu32.lib ${SIMULATION_LINK_LIBRARIES})
  set(SIMULATION_DEPENDENCIES freeglut ${SIMULATION_DEPENDENCIES})
else()
  find_package(GLUT REQUIRED)
  find_package(OpenGL REQUIRED)

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
include_directories(${PROJECT_PATH}/extern/install/NeighborhoodSearch/include)
set(SIMULATION_DEPENDENCIES ${SIMULATION_DEPENDENCIES} Ext_NeighborhoodSearch)
set(SIMULATION_LINK_LIBRARIES ${SIMULATION_LINK_LIBRARIES}
  ${NEIGBORHOOD_SEARCH_LINK_DEPENDENCIES}
	optimized ${NeighborhoodAssemblyName}
	debug ${NeighborhoodAssemblyName}_d)
link_directories(${PROJECT_PATH}/extern/install/NeighborhoodSearch/lib)

############################################################
# GenericParameters
############################################################
include_directories(${PROJECT_PATH}/extern/install/GenericParameters/include)
set(SIMULATION_DEPENDENCIES ${SIMULATION_DEPENDENCIES} Ext_GenericParameters)


add_executable(StaticBoundarySimulator
	main.cpp

	${PROJECT_PATH}/Simulators/Common/SimulatorBase.cpp
	${PROJECT_PATH}/Simulators/Common/SimulatorBase.h
	${PROJECT_PATH}/Simulators/Common/TweakBarParameters.cpp
	${PROJECT_PATH}/Simulators/Common/TweakBarParameters.h

	${VIS_FILES}
)

if(DL_OUTPUT)
	add_definitions(-DDL_OUTPUT)
endif()

add_definitions(-DTW_NO_LIB_PRAGMA -DTW_STATIC)

include_directories(${PROJECT_PATH}/extern/freeglut/include)
include_directories(${PROJECT_PATH}/extern/glew/include)

set_target_properties(StaticBoundarySimulator PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(StaticBoundarySimulator PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(StaticBoundarySimulator PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
add_dependencies(StaticBoundarySimulator ${SIMULATION_DEPENDENCIES})
target_link_libraries(StaticBoundarySimulator ${SIMULATION_LINK_LIBRARIES})
VIS_SOURCE_GROUPS()

set_target_properties(StaticBoundarySimulator PROPERTIES FOLDER "Simulators")

add_custom_command(TARGET StaticBoundarySimulator PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                       ${PROJECT_PATH}/data/shaders $<TARGET_FILE_DIR:StaticBoundarySimulator>/resources/shaders)
