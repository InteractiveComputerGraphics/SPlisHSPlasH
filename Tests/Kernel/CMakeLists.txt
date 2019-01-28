find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )
include_directories(${PROJECT_PATH}/extern/Catch2)

add_executable(KernelTests
	  KernelTests.cpp
)


set_target_properties(KernelTests PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
set_target_properties(KernelTests PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
set_target_properties(KernelTests PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
add_dependencies(KernelTests SPlisHSPlasH Utilities)
target_link_libraries(KernelTests SPlisHSPlasH Utilities)

set_target_properties(KernelTests PROPERTIES FOLDER "Tests")

