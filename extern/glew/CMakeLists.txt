add_library(glew SHARED
	src/glew.c
	include/GL/glew.h
	include/GL/glxew.h
	include/GL/wglew.h
)
target_compile_definitions(glew PRIVATE -DGLEW_NO_GLU -DGLEW_BUILD)
target_include_directories(glew PUBLIC include)

if (WIN32)
	target_link_libraries(glew PUBLIC opengl32.lib)
elseif(APPLE)
   find_library(COCOA_LIBRARY Cocoa)
   find_library(GLUT_LIBRARY GLUT )
   find_library(OpenGL_LIBRARY OpenGL )
   mark_as_advanced(COCOA_LIBRARY GLUT_LIBRARY OpenGL_LIBRARY)
	 target_include_directories(glew PUBLIC /System/Library/Frameworks )
   target_link_libraries(glew PUBLIC ${COCOA_LIBRARY} ${GLUT_LIBRARY} ${OpenGL_LIBRARY})
endif()

set_target_properties(glew PROPERTIES FOLDER "External Dependencies")
set_target_properties(glew PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(glew PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(glew PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
