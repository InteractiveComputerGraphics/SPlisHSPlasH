if (NOT SPH_LIBS_ONLY)
	include_directories(include)

	add_library(glew SHARED
		src/glew.c
		include/GL/glew.h
		include/GL/glxew.h
		include/GL/wglew.h
	)

	set_target_properties(glew PROPERTIES FOLDER "External Dependencies")
	add_definitions(-DGLEW_NO_GLU -DGLEW_BUILD)

	if (WIN32)
		target_link_libraries(glew opengl32.lib)                     
	endif()

	if(APPLE)   
	   include_directories ( /System/Library/Frameworks )
	   find_library(COCOA_LIBRARY Cocoa)
	   find_library(GLUT_LIBRARY GLUT )
	   find_library(OpenGL_LIBRARY OpenGL )
	   mark_as_advanced (COCOA_LIBRARY
						 GLUT_LIBRARY
						 OpenGL_LIBRARY)
	   set(EXTRA_LIBS ${COCOA_LIBRARY} ${GLUT_LIBRARY} ${OpenGL_LIBRARY})
	   target_link_libraries(glew ${EXTRA_LIBS})
	endif (APPLE)

	set_target_properties(glew PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
	set_target_properties(glew PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
	set_target_properties(glew PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
endif()	