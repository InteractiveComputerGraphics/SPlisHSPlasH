add_executable(RigidBodyTests
	  RigidBodyTests.cpp
)
add_dependencies(RigidBodyTests SimulatorBase)
target_include_directories(RigidBodyTests PRIVATE ${PROJECT_SOURCE_DIR}/extern/Catch2)
target_link_libraries(RigidBodyTests SimulatorBase)

set_target_properties(RigidBodyTests PROPERTIES FOLDER "Tests")
set_target_properties(RigidBodyTests PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(RigidBodyTests PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(RigidBodyTests PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
