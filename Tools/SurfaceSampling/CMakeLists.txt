add_executable(SurfaceSampling
	main.cpp
)
add_dependencies(SurfaceSampling SPlisHSPlasH)
target_link_libraries(SurfaceSampling PRIVATE SPlisHSPlasH)

set_target_properties(SurfaceSampling PROPERTIES FOLDER "Tools")
set_target_properties(SurfaceSampling PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(SurfaceSampling PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(SurfaceSampling PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
