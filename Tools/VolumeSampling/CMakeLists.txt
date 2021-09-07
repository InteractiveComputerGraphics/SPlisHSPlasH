add_executable(VolumeSampling
	main.cpp
	SamplingBase.h
	SamplingBase.cpp
	SPHSamplingBase.cpp
	SPHSamplingBase.h
	SPHVolumeSampling.cpp
	SPHVolumeSampling.h
	SPHVolumeSampling_Jiang2015.cpp
	SPHVolumeSampling_Jiang2015.h
)
add_dependencies(VolumeSampling SPlisHSPlasH)
target_link_libraries(VolumeSampling PRIVATE SPlisHSPlasH)

if(DL_OUTPUT)
	target_compile_definitions(VolumeSampling PUBLIC -DDL_OUTPUT)
endif()

set_target_properties(VolumeSampling PROPERTIES FOLDER "Tools")
set_target_properties(VolumeSampling PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(VolumeSampling PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(VolumeSampling PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})

	