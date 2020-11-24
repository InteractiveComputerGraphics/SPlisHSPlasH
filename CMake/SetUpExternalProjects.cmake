if(BUILD_SHARED_LIBS)
	set(LIB_PREFIX ${CMAKE_SHARED_LIBRARY_PREFIX})
	set(LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else(BUILD_SHARED_LIBS)
	set(LIB_PREFIX ${CMAKE_STATIC_LIBRARY_PREFIX})
	set(LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif(BUILD_SHARED_LIBS)

## NeighborhoodSearch
include(NeighborhoodSearch)

## Discregrid
ExternalProject_Add(
	Ext_Discregrid
	PREFIX "${CMAKE_BINARY_DIR}/extern/Discregrid"
	GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/Discregrid.git
	GIT_TAG "0b69062ff9c56fbb6dcecd296652028bedbacf0e"
	INSTALL_DIR ${ExternalInstallDir}/Discregrid
	CMAKE_ARGS -DCMAKE_BUILD_TYPE:STRING=${EXT_CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/Discregrid -DBUILD_CMD_EXECUTABLE:BOOL=0 -DEIGEN3_INCLUDE_DIR:PATH=${EIGEN3_INCLUDE_DIR} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
)
ExternalProject_Get_Property(Ext_Discregrid INSTALL_DIR)
set(DISCREGRID_INCLUDE_DIR ${INSTALL_DIR}/include)
set(DISCREGRID_DEBUG_LIB ${INSTALL_DIR}/lib/${LIB_PREFIX}Discregrid_d${LIB_SUFFIX})
set(DISCREGRID_LIB ${INSTALL_DIR}/lib/${LIB_PREFIX}Discregrid${LIB_SUFFIX})
set(DISCREGRID_LIBRARIES
	optimized ${DISCREGRID_LIB}
	debug ${DISCREGRID_DEBUG_LIB}
)
unset(INSTALL_DIR)

## GenericParameters
ExternalProject_Add(
	Ext_GenericParameters
	PREFIX "${CMAKE_BINARY_DIR}/extern/GenericParameters"
	GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/GenericParameters.git
	GIT_TAG "89ae733fb8b8b9df25d0e44a0e49e51136901e8c"
	INSTALL_DIR ${ExternalInstallDir}/GenericParameters
	CMAKE_ARGS -DCMAKE_BUILD_TYPE=${EXT_CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/GenericParameters -DGENERICPARAMETERS_NO_TESTS:BOOL=1 -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
)
ExternalProject_Get_Property(Ext_GenericParameters INSTALL_DIR)
set(GENERICPARAMETERS_INCLUDE_DIR ${INSTALL_DIR}/include)
unset(INSTALL_DIR)

## PositionBasedDynamics
include(ExternalProject)
ExternalProject_Add(
	Ext_PBD
	PREFIX "${CMAKE_BINARY_DIR}/extern/PositionBasedDynamics"
	GIT_REPOSITORY https://github.com/InteractiveComputerGraphics/PositionBasedDynamics.git
	GIT_TAG "67cea4478d58de55f64f78fcfea16629bfb79152"
	INSTALL_DIR ${ExternalInstallDir}/PositionBasedDynamics
	DEPENDS Ext_GenericParameters Ext_Discregrid
	CMAKE_ARGS -DCMAKE_BUILD_TYPE=${EXT_CMAKE_BUILD_TYPE}
	-DCMAKE_INSTALL_PREFIX:PATH=${ExternalInstallDir}/PositionBasedDynamics
	-DPBD_NO_DEMOS:BOOL=1
	-DPBD_EXTERNALINSTALLDIR:PATH=${ExternalInstallDir}
	-DUSE_DOUBLE_PRECISION:BOOL=${USE_DOUBLE_PRECISION}
	-DGenericParameters_INCLUDE_DIR:PATH=${GENERICPARAMETERS_INCLUDE_DIR}
	-DDiscregrid_INCLUDE_DIR:PATH=${DISCREGRID_INCLUDE_DIR}
	-DDiscregrid_DEBUG_LIB:FILEPATH=${DISCREGRID_DEBUG_LIB}
	-DDiscregrid_LIB:FILEPATH=${DISCREGRID_LIB}
	-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
)
ExternalProject_Get_Property(Ext_PBD INSTALL_DIR)
set(PBD_INCLUDE_DIR ${INSTALL_DIR}/include)
set(PBD_LIBRARIES
	optimized ${INSTALL_DIR}/lib/${LIB_PREFIX}Simulation${LIB_SUFFIX}
	debug ${INSTALL_DIR}/lib/${LIB_PREFIX}Simulation_d${LIB_SUFFIX}
	optimized	${INSTALL_DIR}/lib/${LIB_PREFIX}PositionBasedDynamics${LIB_SUFFIX}
	debug ${INSTALL_DIR}/lib/${LIB_PREFIX}PositionBasedDynamics_d${LIB_SUFFIX}
	optimized ${INSTALL_DIR}/lib/${LIB_PREFIX}Utils${LIB_SUFFIX}
	debug ${INSTALL_DIR}/lib/${LIB_PREFIX}Utils_d${LIB_SUFFIX}
)
unset(INSTALL_DIR)
# This is a workaround to prepend the include directories of PBD.
# TODO change include path of PBD
add_library(PBD_includes INTERFACE)
add_dependencies(PBD_includes Ext_PBD)
target_include_directories(PBD_includes INTERFACE ${PBD_INCLUDE_DIR})
