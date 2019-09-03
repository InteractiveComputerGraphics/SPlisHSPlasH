add_library(tinyexpr STATIC
	tinyexpr.c
	tinyexpr.h
	)

set_target_properties(tinyexpr PROPERTIES FOLDER "External Dependencies")
set_target_properties(tinyexpr PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(tinyexpr PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(tinyexpr PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
