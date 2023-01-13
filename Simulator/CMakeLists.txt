include(${PROJECT_SOURCE_DIR}/GUI/OpenGL/CMakeLists.txt)

include(${PROJECT_SOURCE_DIR}/GUI/imgui/CMakeLists.txt)
set(SIMULATION_LINK_LIBRARIES imgui ${NFD_LIB} SPlisHSPlasH ${PBD_LIBRARIES} PBD_includes)
set(SIMULATION_DEPENDENCIES imgui ${NFD_LIB} SPlisHSPlasH Ext_PBD CopySPlisHSPlasHShaders CopyPBDShaders CopyEmitterModels CopyImguiFonts)
	
set(GUI_SOURCE_FILES 
	${PROJECT_SOURCE_DIR}/Simulator/GUI/imgui/Simulator_GUI_imgui.cpp
	${PROJECT_SOURCE_DIR}/Simulator/GUI/imgui/LogWindow.cpp
	PositionBasedDynamicsWrapper/PBD_Simulator_GUI_imgui.cpp
)
set(GUI_HEADER_FILES
	${PROJECT_SOURCE_DIR}/Simulator/GUI/imgui/Simulator_GUI_imgui.h
	${PROJECT_SOURCE_DIR}/Simulator/GUI/imgui/LogWindow.h
	PositionBasedDynamicsWrapper/PBD_Simulator_GUI_imgui.h
)
	
source_group("Header Files\\GUI\\imgui" FILES ${GUI_HEADER_FILES})
source_group("Source Files\\GUI\\imgui" FILES ${GUI_SOURCE_FILES})


if(WIN32)
  set(SIMULATION_LINK_LIBRARIES opengl32.lib glu32.lib ${SIMULATION_LINK_LIBRARIES})
  set(SIMULATION_DEPENDENCIES ${SIMULATION_DEPENDENCIES})
else()
  set(OpenGL_GL_PREFERENCE GLVND)
  find_package(OpenGL REQUIRED)

  include_directories( ${OPENGL_INCLUDE_DIRS})

  set(SIMULATION_LINK_LIBRARIES
	${SIMULATION_LINK_LIBRARIES}
	${OPENGL_LIBRARIES}
  )
endif()

if (USE_EMBEDDED_PYTHON)
	set(SIMULATION_LINK_LIBRARIES pysplishsplash_embed ${SIMULATION_LINK_LIBRARIES})
	set(SIMULATION_DEPENDENCIES pysplishsplash_embed ${SIMULATION_DEPENDENCIES})
endif()

set(PBDWRAPPER_SOURCE_FILES	
	PositionBasedDynamicsWrapper/PBDWrapper.cpp
	PositionBasedDynamicsWrapper/PBDBoundarySimulator.cpp
)
set(PBDWRAPPER_HEADER_FILES
	PositionBasedDynamicsWrapper/PBDWrapper.h
	PositionBasedDynamicsWrapper/PBDRigidBody.h
	PositionBasedDynamicsWrapper/PBDBoundarySimulator.h
)

set(EXPORTER_SOURCE_FILES
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/ParticleExporter_Partio.cpp
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/ParticleExporter_VTK.cpp
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/RigidBodyExporter_BIN.cpp
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/RigidBodyExporter_OBJ.cpp
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/RigidBodyExporter_PLY.cpp
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/RigidBodyExporter_VTK.cpp
)

set(EXPORTER_HEADER_FILES
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/ExporterBase.h
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/ParticleExporter_Partio.h
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/ParticleExporter_VTK.h
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/RigidBodyExporter_BIN.h
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/RigidBodyExporter_OBJ.h
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/RigidBodyExporter_PLY.h
	${PROJECT_SOURCE_DIR}/Simulator/Exporter/RigidBodyExporter_VTK.h
)

if (USE_EMBEDDED_PYTHON)
	set(EMBEDDED_PYTHON_SOURCE_FILES
		${PROJECT_SOURCE_DIR}/Simulator/ScriptObject.cpp
	)
	set(EMBEDDED_PYTHON_HEADER_FILES
		${PROJECT_SOURCE_DIR}/Simulator/ScriptObject.h
	)
endif(USE_EMBEDDED_PYTHON)

add_library(SimulatorBase
	${PROJECT_SOURCE_DIR}/Simulator/ExporterRegistration.cpp
	${PROJECT_SOURCE_DIR}/Simulator/SceneConfiguration.cpp
	${PROJECT_SOURCE_DIR}/Simulator/SceneConfiguration.h
	${PROJECT_SOURCE_DIR}/Simulator/SimulatorBase.cpp
	${PROJECT_SOURCE_DIR}/Simulator/SimulatorBase.h
	${PROJECT_SOURCE_DIR}/Simulator/BoundarySimulator.cpp
	${PROJECT_SOURCE_DIR}/Simulator/BoundarySimulator.h
	${PROJECT_SOURCE_DIR}/Simulator/StaticBoundarySimulator.cpp
	${PROJECT_SOURCE_DIR}/Simulator/StaticBoundarySimulator.h
	
	${PROJECT_SOURCE_DIR}/Simulator/GUI/OpenGL/Simulator_OpenGL.cpp
	${PROJECT_SOURCE_DIR}/Simulator/GUI/OpenGL/Simulator_OpenGL.h

	${PROJECT_SOURCE_DIR}/Simulator/GUI/Simulator_GUI_Base.h
	${GUI_SOURCE_FILES}
	${GUI_HEADER_FILES}

	${VIS_FILES}
	${GUI_FILES}
	
	${EXPORTER_SOURCE_FILES}
	${EXPORTER_HEADER_FILES}

	${PBDWRAPPER_SOURCE_FILES}
	${PBDWRAPPER_HEADER_FILES}
	
	${EMBEDDED_PYTHON_SOURCE_FILES}
	${EMBEDDED_PYTHON_HEADER_FILES}
)

if(DL_OUTPUT)
	target_compile_definitions(SimulatorBase PUBLIC -DDL_OUTPUT)
endif()


set_target_properties(SimulatorBase PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(SimulatorBase PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(SimulatorBase PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
add_dependencies(SimulatorBase ${SIMULATION_DEPENDENCIES})
target_link_libraries(SimulatorBase PUBLIC ${SIMULATION_LINK_LIBRARIES})
VIS_SOURCE_GROUPS()
GUI_SOURCE_GROUPS()

source_group("Header Files\\Exporter" FILES ${EXPORTER_HEADER_FILES})
source_group("Source Files\\Exporter" FILES ${EXPORTER_SOURCE_FILES})
source_group("Header Files\\PBD" FILES ${PBDWRAPPER_HEADER_FILES})
source_group("Source Files\\PBD" FILES ${PBDWRAPPER_SOURCE_FILES})
source_group("Header Files\\GUI\\OpenGL" FILES ${PROJECT_SOURCE_DIR}/Simulator/GUI/OpenGL/Simulator_OpenGL.h)
source_group("Source Files\\GUI\\OpenGL" FILES ${PROJECT_SOURCE_DIR}/Simulator/GUI/OpenGL/Simulator_OpenGL.cpp)

add_executable(SPHSimulator main.cpp)
target_link_libraries(SPHSimulator SimulatorBase)
add_dependencies(SPHSimulator SimulatorBase)

set_target_properties(SPHSimulator PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(SPHSimulator PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(SPHSimulator PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})

