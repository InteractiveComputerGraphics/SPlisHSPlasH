find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

add_executable(SurfaceSampling
	main.cpp

	CMakeLists.txt
)

set_target_properties(SurfaceSampling PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(SurfaceSampling PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(SurfaceSampling PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
add_dependencies(SurfaceSampling SPlisHSPlasH Utilities partio zlib)
target_link_libraries(SurfaceSampling SPlisHSPlasH Utilities partio zlib)

source_group("Header Files\\Utils" FILES ${UTILS_HEADER_FILES})
source_group("Source Files\\Utils" FILES ${UTILS_SOURCE_FILES})

set_target_properties(SurfaceSampling PROPERTIES FOLDER "Tools")
