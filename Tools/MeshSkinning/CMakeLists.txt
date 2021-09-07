add_executable(MeshSkinning
	main.cpp
)
add_dependencies(MeshSkinning SPlisHSPlasH)
target_link_libraries(MeshSkinning PRIVATE SPlisHSPlasH)

set_target_properties(MeshSkinning PROPERTIES FOLDER "Tools")
set_target_properties(MeshSkinning PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(MeshSkinning PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(MeshSkinning PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
