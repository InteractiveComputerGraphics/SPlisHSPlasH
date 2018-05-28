find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

add_executable(SurfaceSampling
	main.cpp
)

set_target_properties(SurfaceSampling PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(SurfaceSampling PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(SurfaceSampling PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
add_dependencies(SurfaceSampling SPlisHSPlasH Utilities partio zlib MD5)
target_link_libraries(SurfaceSampling SPlisHSPlasH Utilities partio zlib MD5)

set_target_properties(SurfaceSampling PROPERTIES FOLDER "Tools")
