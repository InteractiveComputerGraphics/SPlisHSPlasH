find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

add_library(Utilities
	
	FileSystem.cpp
	FileSystem.h
	OBJLoader.cpp
	OBJLoader.h
	PartioReaderWriter.cpp
	PartioReaderWriter.h
	SceneLoader.cpp
	SceneLoader.h
	StringTools.cpp
	StringTools.h
	
	CMakeLists.txt
)

add_dependencies(Utilities partio)
