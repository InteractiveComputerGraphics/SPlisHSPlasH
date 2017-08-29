include(${PROJECT_PATH}/Visualization/CMakeLists.txt)
add_definitions(-DPBD_DATA_PATH="../data")

subdirs(PartioViewer SurfaceSampling VolumeSampling)

