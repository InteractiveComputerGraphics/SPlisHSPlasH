set(IMGUI_HEADER_FILES
	${PROJECT_SOURCE_DIR}/GUI/imgui/imguiParameters.h
	) 

	
set(IMGUI_SOURCE_FILES
	${PROJECT_SOURCE_DIR}/GUI/imgui/imguiParameters.cpp   
	) 
	

set(GUI_FILES
	${IMGUI_HEADER_FILES}
	${IMGUI_SOURCE_FILES}
	) 
	
# A macro for the source groups
macro(GUI_SOURCE_GROUPS)
	source_group("Header Files\\GUI\\imgui" FILES ${IMGUI_HEADER_FILES})
	source_group("Source Files\\GUI\\imgui" FILES ${IMGUI_SOURCE_FILES})
endmacro(GUI_SOURCE_GROUPS)



