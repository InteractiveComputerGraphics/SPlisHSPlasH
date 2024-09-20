include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_LOCAL_CHANGES)
if (${GIT_LOCAL_CHANGES} STREQUAL "DIRTY")
	set(COMPILER_MESSAGE "#pragma WARNING(Local changes not committed.)")
endif()

file (STRINGS "../version.txt" SPLISHSPLASH_VERSION)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/Version.h @ONLY)




add_library(Utilities
	
	BinaryFileReaderWriter.h
	Counting.h
	FileSystem.h
	Logger.h
	OBJLoader.h
	PartioReaderWriter.cpp
	PartioReaderWriter.h
	PLYLoader.h
	StringTools.h	
	SystemInfo.h
	Timing.h
	Version.h
)

add_dependencies(Utilities partio zlib)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	target_compile_options(Utilities PUBLIC -idirafter${TOPLEVEL_INCLUDE_DIR} -idirafter${EIGEN3_INCLUDE_DIR})
else()
	target_include_directories(Utilities PUBLIC ${TOPLEVEL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR})
endif()
target_link_libraries(Utilities INTERFACE partio zlib MD5)
