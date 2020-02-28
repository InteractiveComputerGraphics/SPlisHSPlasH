set(GUI_HEADER_FILES
	${PROJECT_SOURCE_DIR}/GUI/TweakBar/TweakBarParameters.h
	) 

	
set(GUI_SOURCE_FILES
	${PROJECT_SOURCE_DIR}/GUI/TweakBar/TweakBarParameters.cpp   
	) 
	

set(GUI_FILES
	${GUI_HEADER_FILES}
	${GUI_SOURCE_FILES}
	) 
	
# A macro for the source groups
macro(GUI_SOURCE_GROUPS)
	source_group("Header Files\\GUI\\TweakBar" FILES ${GUI_HEADER_FILES})
	source_group("Source Files\\GUI\\TweakBar" FILES ${GUI_SOURCE_FILES})
endmacro(GUI_SOURCE_GROUPS)



