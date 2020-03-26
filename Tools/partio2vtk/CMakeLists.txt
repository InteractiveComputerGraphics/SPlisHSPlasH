add_executable(partio2vtk
	main.cpp
)
add_dependencies(partio2vtk SPlisHSPlasH)
target_link_libraries(partio2vtk PRIVATE SPlisHSPlasH)

set_target_properties(partio2vtk PROPERTIES FOLDER "Tools")
set_target_properties(partio2vtk PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(partio2vtk PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(partio2vtk PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
