if (NOT SPH_LIBS_ONLY)
	include(${PROJECT_PATH}/Visualization/CMakeLists.txt)
	add_definitions(-DPBD_DATA_PATH="../data")

	add_subdirectory(PartioViewer)
	add_subdirectory(SurfaceSampling)
	add_subdirectory(VolumeSampling)
	add_subdirectory(partio2vtk)
endif()	

