set(SIMULATION_LINK_LIBRARIES Utilities AntTweakBar glew partio zlib SPlisHSPlasH MD5)
set(SIMULATION_DEPENDENCIES Utilities AntTweakBar glew partio zlib SPlisHSPlasH MD5)

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
# CompactNSearch
############################################################
include_directories(${PROJECT_PATH}/extern/install/CompactNSearch/include)
set(SIMULATION_LINK_LIBRARIES ${SIMULATION_LINK_LIBRARIES} 
	optimized CompactNSearch 
	debug CompactNSearch_d)
set(SIMULATION_DEPENDENCIES ${SIMULATION_DEPENDENCIES} ExternalProject_CompactNSearch)
link_directories(${PROJECT_PATH}/extern/install/CompactNSearch/lib)

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

set(SIMULATION_DEPENDENCIES ExternalProject_PositionBasedDynamics ${SIMULATION_DEPENDENCIES})
link_directories(${PROJECT_PATH}/extern/install/PositionBasedDynamics/lib)

add_executable(DynamicBoundaryDemo
	main.cpp
	
	${PROJECT_PATH}/Demos/Common/DemoBase.cpp
	${PROJECT_PATH}/Demos/Common/DemoBase.h
  
	${VIS_FILES}          

	${PBDWRAPPER_SOURCE_FILES}
	${PBDWRAPPER_HEADER_FILES}
  
	CMakeLists.txt
)

add_definitions(-DTW_NO_LIB_PRAGMA -DTW_STATIC)

include_directories(${PROJECT_PATH}/extern/freeglut/include)
include_directories(${PROJECT_PATH}/extern/glew/include)

set_target_properties(DynamicBoundaryDemo PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(DynamicBoundaryDemo PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(DynamicBoundaryDemo PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
add_dependencies(DynamicBoundaryDemo ${SIMULATION_DEPENDENCIES})
target_link_libraries(DynamicBoundaryDemo ${SIMULATION_LINK_LIBRARIES})
VIS_SOURCE_GROUPS()

source_group("Header Files\\PBD" FILES ${PBDWRAPPER_HEADER_FILES})
source_group("Source Files\\PBD" FILES ${PBDWRAPPER_SOURCE_FILES})

set_target_properties(DynamicBoundaryDemo PROPERTIES FOLDER "Demos")


