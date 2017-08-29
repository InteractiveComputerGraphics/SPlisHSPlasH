find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

add_executable(VolumeSampling
	main.cpp
	WindingNumbers.cpp
	WindingNumbers.h
	
	CMakeLists.txt
)

set_target_properties(VolumeSampling PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(VolumeSampling PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(VolumeSampling PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
add_dependencies(VolumeSampling SPlisHSPlasH Utilities partio zlib MD5)
target_link_libraries(VolumeSampling SPlisHSPlasH Utilities partio zlib MD5)

set_target_properties(VolumeSampling PROPERTIES FOLDER "Tools")