set(VIS_HEADER_FILES
	${PROJECT_SOURCE_DIR}/GUI/OpenGL/MiniGL.h
	${PROJECT_SOURCE_DIR}/GUI/OpenGL/Selection.h
	${PROJECT_SOURCE_DIR}/GUI/OpenGL/Shader.h
	${PROJECT_SOURCE_DIR}/GUI/OpenGL/colormaps/colormap_jet.h
	${PROJECT_SOURCE_DIR}/GUI/OpenGL/colormaps/colormap_plasma.h
	) 

	
set(VIS_SOURCE_FILES
	${PROJECT_SOURCE_DIR}/extern/glfw/deps/glad_gl.c
	${PROJECT_SOURCE_DIR}/GUI/OpenGL/MiniGL.cpp
	${PROJECT_SOURCE_DIR}/GUI/OpenGL/Shader.cpp
	) 
	

set(VIS_FILES
	${VIS_HEADER_FILES}
	${VIS_SOURCE_FILES}
	) 
	
# A macro for the source groups
macro(VIS_SOURCE_GROUPS)
	source_group("Header Files\\GUI\\OpenGL" FILES ${VIS_HEADER_FILES})
	source_group("Source Files\\GUI\\OpenGL" FILES ${VIS_SOURCE_FILES})
endmacro(VIS_SOURCE_GROUPS)



