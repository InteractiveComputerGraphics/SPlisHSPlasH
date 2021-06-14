add_executable(FoamGenerator
	main.cpp
	
	FoamKernel.cpp
	FoamKernel.h
)

add_dependencies(FoamGenerator SPlisHSPlasH)
target_link_libraries(FoamGenerator PRIVATE SPlisHSPlasH)

set_target_properties(FoamGenerator PROPERTIES FOLDER "Tools")
set_target_properties(FoamGenerator PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(FoamGenerator PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(FoamGenerator PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})

