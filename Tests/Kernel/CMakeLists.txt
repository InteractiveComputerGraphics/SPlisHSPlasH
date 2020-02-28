add_executable(KernelTests
	  KernelTests.cpp
)
add_dependencies(KernelTests SPlisHSPlasH)
target_include_directories(KernelTests PRIVATE ${PROJECT_SOURCE_DIR}/extern/Catch2)
target_link_libraries(KernelTests PRIVATE SPlisHSPlasH)

set_target_properties(KernelTests PROPERTIES FOLDER "Tests")
set_target_properties(KernelTests PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(KernelTests PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(KernelTests PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
