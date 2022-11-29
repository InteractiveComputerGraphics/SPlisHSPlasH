include(${PROJECT_SOURCE_DIR}/GUI/OpenGL/CMakeLists.txt)

include(${PROJECT_SOURCE_DIR}/GUI/imgui/CMakeLists.txt)
set(VIEWER_LINK_LIBRARIES imgui SPlisHSPlasH)
set(VIEWER_DEPENDENCIES imgui SPlisHSPlasH CopySPlisHSPlasHShaders CopyImguiFonts)

set(GUI_SOURCE_FILES 
	GUI/imgui/PartioViewer_GUI_imgui.cpp
)
set(GUI_HEADER_FILES
	GUI/imgui/PartioViewer_GUI_imgui.h
)

source_group("Header Files\\GUI\\imgui" FILES ${GUI_HEADER_FILES})
source_group("Source Files\\GUI\\imgui" FILES ${GUI_SOURCE_FILES})

if(WIN32)
  set(VIEWER_LINK_LIBRARIES opengl32.lib glu32.lib ${VIEWER_LINK_LIBRARIES})
  set(VIEWER_DEPENDENCIES ${VIEWER_DEPENDENCIES})
else()
  set(OpenGL_GL_PREFERENCE GLVND)
  find_package(OpenGL REQUIRED)

  include_directories( ${OPENGL_INCLUDE_DIRS})

  set(VIEWER_LINK_LIBRARIES
    ${VIEWER_LINK_LIBRARIES}
    ${OPENGL_LIBRARIES}
  )
endif()



add_executable(PartioViewer
	main.cpp
	
	PartioViewer.cpp
	PartioViewer.h
	
	GUI/OpenGL/PartioViewer_OpenGL.cpp
	GUI/OpenGL/PartioViewer_OpenGL.h
	
	GUI/PartioViewer_GUI_Base.h
	${GUI_SOURCE_FILES}
	${GUI_HEADER_FILES}
	
	${PROJECT_SOURCE_DIR}/extern/toojpeg/toojpeg.cpp
	${PROJECT_SOURCE_DIR}/extern/toojpeg/toojpeg.h
	
	${VIS_FILES}
	${GUI_FILES}
)

add_dependencies(PartioViewer ${VIEWER_DEPENDENCIES})
target_link_libraries(PartioViewer ${VIEWER_LINK_LIBRARIES})
VIS_SOURCE_GROUPS()
GUI_SOURCE_GROUPS()

set_target_properties(PartioViewer PROPERTIES FOLDER "Tools")
set_target_properties(PartioViewer PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(PartioViewer PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(PartioViewer PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})

source_group("Header Files\\GUI" FILES 	GUI/PartioViewer_GUI_Base.h)
source_group("Header Files\\GUI\\OpenGL" FILES 	GUI/OpenGL/PartioViewer_OpenGL.h)
source_group("Source Files\\GUI\\OpenGL" FILES 	GUI/OpenGL/PartioViewer_OpenGL.cpp)
