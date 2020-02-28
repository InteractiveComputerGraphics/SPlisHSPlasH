add_executable(VolumeSampling
	main.cpp
)
add_dependencies(VolumeSampling SPlisHSPlasH)
target_link_libraries(VolumeSampling PRIVATE SPlisHSPlasH)

set_target_properties(VolumeSampling PROPERTIES FOLDER "Tools")
set_target_properties(VolumeSampling PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(VolumeSampling PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(VolumeSampling PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
