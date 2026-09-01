include(ExternalProject)

set(NeighborhoodSearch "CompactNSearch" CACHE STRING "NeighborhoodSearch chosen by the user at CMake configure time")
set_property(CACHE NeighborhoodSearch PROPERTY STRINGS CompactNSearch cuNSearch TreeNSearch)


if ("${NeighborhoodSearch}" STREQUAL "cuNSearch")

	message(STATUS "Use neighborhood search: cuNSearch")
	
	add_definitions( -DUSE_cuNSearch)	

	if(USE_DOUBLE_PRECISION)
		message("Use cuNSearch with single precision to get a better performance.")
	endif()

	enable_language(CUDA)
	find_package(CUDA 9.0 REQUIRED)
	## cuNSearch
	ExternalProject_Add(
	   Ext_NeighborhoodSearch
	   PREFIX "${CMAKE_SOURCE_DIR}/extern/cuNSearch"
	   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/cuNSearch.git
	   GIT_TAG "aba3da18cb4f45cd05d729465d1725891ffc33da"
	   INSTALL_DIR ${ExternalInstallDir}/NeighborhoodSearch
	   CMAKE_ARGS -DCMAKE_BUILD_TYPE=${EXT_CMAKE_BUILD_TYPE} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE} -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/NeighborhoodSearch -DCUNSEARCH_USE_DOUBLE_PRECISION:BOOL=${USE_DOUBLE_PRECISION} -DBUILD_DEMO:BOOL=OFF
	   )

	set(NEIGHBORHOOD_ASSEMBLY_NAME cuNSearch)
	set(NEIGBORHOOD_SEARCH_LINK_DEPENDENCIES general ${CUDA_LIBRARIES})
	add_compile_options(-DGPU_NEIGHBORHOOD_SEARCH)

elseif ("${NeighborhoodSearch}" STREQUAL "CompactNSearch")
	
	message(STATUS "Use neighborhood search: CompactNSearch")
	
	add_definitions( -DUSE_CompactNSearch)	

	## CompactNSearch
	ExternalProject_Add(
	   Ext_NeighborhoodSearch
	   PREFIX "${CMAKE_BINARY_DIR}/extern/CompactNSearch"
	   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/CompactNSearch.git
	   GIT_TAG "b40afcf47fe1963b363eba2371f04b42720fcb1d"
	   INSTALL_DIR ${ExternalInstallDir}/NeighborhoodSearch
	   CMAKE_ARGS -DCMAKE_BUILD_TYPE=${EXT_CMAKE_BUILD_TYPE} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE} -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/NeighborhoodSearch -DUSE_DOUBLE_PRECISION:BOOL=${USE_DOUBLE_PRECISION} -DBUILD_DEMO:BOOL=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.10 -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
	)
	set(NEIGHBORHOOD_ASSEMBLY_NAME CompactNSearch)
	if(WIN32)
		add_definitions( -D_SILENCE_STDEXT_ARR_ITERS_DEPRECATION_WARNING)
	endif() 

elseif ("${NeighborhoodSearch}" STREQUAL "TreeNSearch")
	
	message(STATUS "Use neighborhood search: TreeNSearch")
	
	add_definitions( -DUSE_TreeNSearch)	
	
	## TreeNSearch
	ExternalProject_Add(
	   Ext_NeighborhoodSearch
	   PREFIX "${CMAKE_BINARY_DIR}/extern/TreeNSearch"
	   GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/TreeNSearch.git
	   GIT_TAG "61cf111f65c01f1d71fb5023ac72ce6bb0e8d6f9"
	   INSTALL_DIR ${ExternalInstallDir}/NeighborhoodSearch
	   CMAKE_ARGS -DCMAKE_BUILD_TYPE=${EXT_CMAKE_BUILD_TYPE} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} -DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE} -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/NeighborhoodSearch -DUSE_DOUBLE_PRECISION:BOOL=${USE_DOUBLE_PRECISION} -DBUILD_DEMO:BOOL=OFF -DCMAKE_DEBUG_POSTFIX=_d -DCMAKE_RELWITHDEBINFO_POSTFIX=_rd -DCMAKE_MINSIZEREL_POSTFIX=_ms
	)
	set(NEIGHBORHOOD_ASSEMBLY_NAME TreeNSearch)
	
endif()

ExternalProject_Get_Property(
	Ext_NeighborhoodSearch
	INSTALL_DIR
)
set(NEIGHBORHOOD_SEARCH_LIBRARIES
	optimized ${INSTALL_DIR}/lib/${LIB_PREFIX}${NEIGHBORHOOD_ASSEMBLY_NAME}${LIB_SUFFIX}
	debug ${INSTALL_DIR}/lib/${LIB_PREFIX}${NEIGHBORHOOD_ASSEMBLY_NAME}_d${LIB_SUFFIX}
	${NEIGBORHOOD_SEARCH_LINK_DEPENDENCIES}
)
set(NEIGHBORHOOD_SEARCH_INCLUDE_DIR ${INSTALL_DIR}/include/)

unset(INSTALL_DIR)
