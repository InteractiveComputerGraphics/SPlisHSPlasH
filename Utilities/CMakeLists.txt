find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )
include_directories(${PROJECT_PATH}/extern/spdlog/include)


include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_LOCAL_CHANGES)
if (${GIT_LOCAL_CHANGES} STREQUAL "DIRTY")
	set(COMPILER_MESSAGE "#pragma WARNING(Local changes not committed.)")
endif()
configure_file("${PROJECT_PATH}/Utilities/Version.h.in" "${PROJECT_PATH}/Utilities/Version.h" @ONLY)




add_library(Utilities
	
	FileSystem.h
	Logger.h
	OBJLoader.h
	PartioReaderWriter.cpp
	PartioReaderWriter.h
	StringTools.h	
	SystemInfo.h
	Timing.h
	Version.h
	
	CMakeLists.txt
)

add_dependencies(Utilities partio)
