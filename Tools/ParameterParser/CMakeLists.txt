add_executable(ParameterParser
	main.cpp
	
	JsonSchemaGenerator.cpp
	JsonSchemaGenerator.h
	ParameterObjectParser.cpp
	ParameterObjectParser.h
	SceneExampleGenerator.cpp
	SceneExampleGenerator.h
)
add_dependencies(ParameterParser SimulatorBase)
target_link_libraries(ParameterParser PRIVATE SimulatorBase)

set_target_properties(ParameterParser PROPERTIES FOLDER "Tools")
set_target_properties(ParameterParser PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(ParameterParser PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(ParameterParser PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
