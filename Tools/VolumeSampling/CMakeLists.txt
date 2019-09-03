find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

set(SIMULATION_LINK_LIBRARIES SPlisHSPlasH Utilities)
set(SIMULATION_DEPENDENCIES SPlisHSPlasH Utilities)

############################################################
# DiscreGrid
############################################################
include_directories(${PROJECT_PATH}/extern/install/Discregrid/include)
set(SIMULATION_LINK_LIBRARIES ${SIMULATION_LINK_LIBRARIES} 
	optimized Discregrid 
	debug Discregrid_d)
set(SIMULATION_DEPENDENCIES ${SIMULATION_DEPENDENCIES} Ext_Discregrid)	
link_directories(${PROJECT_PATH}/extern/install/Discregrid/lib)


add_executable(VolumeSampling
	main.cpp
)

add_dependencies(VolumeSampling ${SIMULATION_DEPENDENCIES})
target_link_libraries(VolumeSampling ${SIMULATION_LINK_LIBRARIES})


set_target_properties(VolumeSampling PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(VolumeSampling PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(VolumeSampling PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
add_dependencies(VolumeSampling SPlisHSPlasH Utilities partio zlib MD5)
target_link_libraries(VolumeSampling SPlisHSPlasH Utilities partio zlib MD5)

set_target_properties(VolumeSampling PROPERTIES FOLDER "Tools")
