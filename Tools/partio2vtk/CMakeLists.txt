set(TARGET_NAME partio2vtk)

find_package( Eigen3 REQUIRED )

add_executable(${TARGET_NAME}
	main.cpp
)

target_link_libraries(${TARGET_NAME} SPlisHSPlasH Utilities partio zlib MD5)

target_include_directories(${TARGET_NAME} PRIVATE
	${EIGEN3_INCLUDE_DIR}
)

set_target_properties(${TARGET_NAME} PROPERTIES
	DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX}
	RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX}
	MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX}
	FOLDER "Tools"
)

add_dependencies(${TARGET_NAME} SPlisHSPlasH Utilities partio zlib MD5)
