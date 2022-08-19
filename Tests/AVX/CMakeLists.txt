add_executable(AVXTests
	  AVXTests.cpp
)
add_dependencies(AVXTests SPlisHSPlasH)
target_include_directories(AVXTests PRIVATE ${PROJECT_SOURCE_DIR}/extern/Catch2)
target_link_libraries(AVXTests PRIVATE SPlisHSPlasH)

set_target_properties(AVXTests PROPERTIES FOLDER "Tests")
set_target_properties(AVXTests PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(AVXTests PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(AVXTests PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
