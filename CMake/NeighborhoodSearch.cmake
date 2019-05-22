include(ExternalProject)

option(USE_GPU_NEIGHBORHOOD_SEARCH "Use GPU neighborhood search" OFF)

if(USE_GPU_NEIGHBORHOOD_SEARCH)
	message(STATUS "Use cuNSearch for neighborhood search")
	
	if(USE_DOUBLE_PRECISION)
		message("Use cuNSearch with single precision to get a better performance.")
	endif()

	find_package(CUDA 9.0 REQUIRED)
	## cuNSearch
	ExternalProject_Add(
	   Ext_NeighborhoodSearch
	   PREFIX "${CMAKE_SOURCE_DIR}/extern/cuNSearch"
	   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/cuNSearch.git
	   GIT_TAG "d4c4eb63648cb243c8634e6bad34d90e7cba74ef"
	   INSTALL_DIR ${ExternalInstallDir}/NeighborhoodSearch
	   CMAKE_ARGS -DCMAKE_BUILD_TYPE=${EXT_CMAKE_BUILD_TYPE} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE} -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/NeighborhoodSearch -DCUNSEARCH_USE_DOUBLE_PRECISION:BOOL=${USE_DOUBLE_PRECISION} -DBUILD_DEMO:BOOL=OFF -DCUDA_USE_STATIC_CUDA_RUNTIME:BOOL=OFF
	   )

	set(CUDA_USE_STATIC_CUDA_RUNTIME OFF CACHE BOOL "" FORCE)
	set(NeighborhoodAssemblyName "cuNSearch")
	set(NEIGBORHOOD_SEARCH_LINK_DEPENDENCIES ${CUDA_LIBRARIES})
	add_compile_options(-DGPU_NEIGHBORHOOD_SEARCH)

else()
	## CompactNSearch
	ExternalProject_Add(
	   Ext_NeighborhoodSearch
	   PREFIX "${CMAKE_SOURCE_DIR}/extern/CompactNSearch"
	   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/CompactNSearch.git
	   GIT_TAG "daf73f4d321ea3cfbdc88d23f56d8aaa0083479a"
	   INSTALL_DIR ${ExternalInstallDir}/NeighborhoodSearch
	   CMAKE_ARGS -DCMAKE_BUILD_TYPE=${EXT_CMAKE_BUILD_TYPE} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE} -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/NeighborhoodSearch -DUSE_DOUBLE_PRECISION:BOOL=${USE_DOUBLE_PRECISION}
	)

	set(NeighborhoodAssemblyName "CompactNSearch")

endif(USE_GPU_NEIGHBORHOOD_SEARCH)

