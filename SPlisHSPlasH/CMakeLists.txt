set(WCSPH_HEADER_FILES
	WCSPH/SimulationDataWCSPH.h
	WCSPH/TimeStepWCSPH.h
	)
	
set(WCSPH_SOURCE_FILES
	WCSPH/SimulationDataWCSPH.cpp
	WCSPH/TimeStepWCSPH.cpp
	)		
	
set(PCISPH_HEADER_FILES
	PCISPH/SimulationDataPCISPH.h
	PCISPH/TimeStepPCISPH.h
	)
	
set(PCISPH_SOURCE_FILES
	PCISPH/SimulationDataPCISPH.cpp
	PCISPH/TimeStepPCISPH.cpp
	)

set(PBF_HEADER_FILES
	PBF/SimulationDataPBF.h
	PBF/TimeStepPBF.h
	PBF/TimeIntegration.h
	)
	
set(PBF_SOURCE_FILES
	PBF/SimulationDataPBF.cpp
	PBF/TimeStepPBF.cpp
	PBF/TimeIntegration.cpp
	)	
	
set(IISPH_HEADER_FILES
	IISPH/SimulationDataIISPH.h
	IISPH/TimeStepIISPH.h
	)
	
set(IISPH_SOURCE_FILES
	IISPH/SimulationDataIISPH.cpp
	IISPH/TimeStepIISPH.cpp
	)		
	
set(DFSPH_HEADER_FILES
	DFSPH/SimulationDataDFSPH.h
	DFSPH/TimeStepDFSPH.h
	)
	
set(DFSPH_SOURCE_FILES
	DFSPH/SimulationDataDFSPH.cpp
	DFSPH/TimeStepDFSPH.cpp
	)	
	
set(PF_HEADER_FILES
	PF/SimulationDataPF.h
	PF/TimeStepPF.h
	)
	
set(PF_SOURCE_FILES
	PF/SimulationDataPF.cpp
	PF/TimeStepPF.cpp
	)
	
set(ICSPH_HEADER_FILES
	ICSPH/SimulationDataICSPH.h
	ICSPH/TimeStepICSPH.h
	)
	
set(ICSPH_SOURCE_FILES
	ICSPH/SimulationDataICSPH.cpp
	ICSPH/TimeStepICSPH.cpp
	)	
	
set(SURFACETENSION_HEADER_FILES
	SurfaceTension/SurfaceTensionBase.h
	SurfaceTension/SurfaceTension_Becker2007.h
	SurfaceTension/SurfaceTension_Akinci2013.h
	SurfaceTension/SurfaceTension_He2014.h
	)
	
set(SURFACETENSION_SOURCE_FILES
	SurfaceTension/SurfaceTensionBase.cpp
	SurfaceTension/SurfaceTension_Becker2007.cpp
	SurfaceTension/SurfaceTension_Akinci2013.cpp
	SurfaceTension/SurfaceTension_He2014.cpp
	)	

set(VISCOSITY_HEADER_FILES
	Viscosity/ViscosityBase.h
	Viscosity/Viscosity_Standard.h
	Viscosity/Viscosity_Bender2017.h
	Viscosity/Viscosity_Peer2015.h
	Viscosity/Viscosity_Peer2016.h
	Viscosity/Viscosity_Takahashi2015.h
	Viscosity/Viscosity_Weiler2018.h
	)
	
set(VISCOSITY_SOURCE_FILES
	Viscosity/ViscosityBase.cpp
	Viscosity/Viscosity_Standard.cpp
	Viscosity/Viscosity_Bender2017.cpp
	Viscosity/Viscosity_Peer2015.cpp
	Viscosity/Viscosity_Peer2016.cpp
	Viscosity/Viscosity_Takahashi2015.cpp
	Viscosity/Viscosity_Weiler2018.cpp
	)	
	
set(VORTICITY_HEADER_FILES
	Vorticity/VorticityBase.h
	Vorticity/MicropolarModel_Bender2017.h
	Vorticity/VorticityConfinement.h
	)
	
set(VORTICITY_SOURCE_FILES
	Vorticity/VorticityBase.cpp
	Vorticity/MicropolarModel_Bender2017.cpp
	Vorticity/VorticityConfinement.cpp
	)	

set(DRAG_HEADER_FILES
	Drag/DragBase.h
	Drag/DragForce_Gissler2017.h
	Drag/DragForce_Macklin2014.h
	)
	
set(DRAG_SOURCE_FILES
	Drag/DragBase.cpp
	Drag/DragForce_Gissler2017.cpp
	Drag/DragForce_Macklin2014.cpp
	)	
	
set(ELASTICITY_HEADER_FILES
	Elasticity/ElasticityBase.h
	Elasticity/Elasticity_Becker2009.h	
	Elasticity/Elasticity_Peer2018.h
	Elasticity/Elasticity_Kugelstadt2021.h
	)
	
set(ELASTICITY_SOURCE_FILES
	Elasticity/ElasticityBase.cpp
	Elasticity/Elasticity_Becker2009.cpp
	Elasticity/Elasticity_Peer2018.cpp
	Elasticity/Elasticity_Kugelstadt2021.cpp
	)
	
set(UTILS_HEADER_FILES
	Utilities/AVX_math.h
	Utilities/CholeskyAVXSolver.h	
	Utilities/DebugTools.h
	Utilities/MathFunctions.h
	Utilities/GaussQuadrature.h
	Utilities/MatrixFreeSolver.h
	Utilities/MeshImport.h
	Utilities/SurfaceSampling.h
	Utilities/PoissonDiskSampling.h
	Utilities/RegularTriangleSampling.h
	Utilities/RegularSampling2D.h
	Utilities/SceneLoader.h
	Utilities/SceneParameterObjects.h
	Utilities/SceneWriter.h
	Utilities/VolumeSampling.h
	Utilities/SDFFunctions.h
	Utilities/SimpleQuadrature.h
	Utilities/WindingNumbers.h
	)
	
set(UTILS_SOURCE_FILES
	Utilities/CholeskyAVXSolver.cpp
	Utilities/DebugTools.cpp
	Utilities/MathFunctions.cpp
	Utilities/GaussQuadrature.cpp
	Utilities/MeshImport.cpp
	Utilities/PoissonDiskSampling.cpp
	Utilities/RegularTriangleSampling.cpp
	Utilities/RegularSampling2D.cpp
	Utilities/SceneLoader.cpp
	Utilities/SceneParameterObjects.cpp
	Utilities/SceneWriter.cpp
	Utilities/VolumeSampling.cpp
	Utilities/SDFFunctions.cpp
	Utilities/SimpleQuadrature.cpp
	Utilities/WindingNumbers.cpp
	)	
	
if (USE_THIRD_PARTY_METHODS)
	list(APPEND SURFACETENSION_HEADER_FILES 
		SurfaceTension/SurfaceTension_ZorillaRitter2020_haltonVec323.h
		SurfaceTension/SurfaceTension_ZorillaRitter2020.h
	)
	
	list(APPEND SURFACETENSION_SOURCE_FILES
		SurfaceTension/SurfaceTension_ZorillaRitter2020.cpp
	)
endif()
	
add_library(SPlisHSPlasH
	Common.h
	
	NeighborhoodSearch.h

	AnimationField.cpp
	AnimationField.h
	AnimationFieldSystem.h
	AnimationFieldSystem.cpp
	BoundaryModel.cpp
	BoundaryModel.h
	BoundaryModel_Akinci2012.cpp
	BoundaryModel_Akinci2012.h
	BoundaryModel_Bender2019.cpp
	BoundaryModel_Bender2019.h
	BoundaryModel_Koschier2017.cpp
	BoundaryModel_Koschier2017.h
	Emitter.cpp
	Emitter.h
	EmitterSystem.cpp
	EmitterSystem.h
	FluidModel.cpp
	FluidModel.h
	RigidBodyObject.h
	NonPressureForceRegistration.cpp
	SPHKernels.cpp
	SPHKernels.h
	StaticRigidBody.h
	Simulation.cpp
	Simulation.h
	TimeManager.cpp
	TimeManager.h	
	TimeStep.cpp
	TimeStep.h
	TriangleMesh.cpp
	TriangleMesh.h
	NonPressureForceBase.cpp
	NonPressureForceBase.h
	XSPH.cpp
	XSPH.h

	
	${WCSPH_HEADER_FILES}
	${WCSPH_SOURCE_FILES}
	
	${PCISPH_HEADER_FILES}
	${PCISPH_SOURCE_FILES}
	
	${PBF_HEADER_FILES}
	${PBF_SOURCE_FILES}
	
	${IISPH_HEADER_FILES}
	${IISPH_SOURCE_FILES}
	
	${DFSPH_HEADER_FILES}
	${DFSPH_SOURCE_FILES}
	
	${PF_HEADER_FILES}
	${PF_SOURCE_FILES}
	
	${ICSPH_HEADER_FILES}
	${ICSPH_SOURCE_FILES}
	
	${VISCOSITY_HEADER_FILES}
	${VISCOSITY_SOURCE_FILES}
	
	${SURFACETENSION_HEADER_FILES}
	${SURFACETENSION_SOURCE_FILES}
	
	${VORTICITY_HEADER_FILES}
	${VORTICITY_SOURCE_FILES}
	
	${DRAG_HEADER_FILES}
	${DRAG_SOURCE_FILES}

	${ELASTICITY_HEADER_FILES}
	${ELASTICITY_SOURCE_FILES}

	${UTILS_HEADER_FILES}
	${UTILS_SOURCE_FILES}
)

add_dependencies(SPlisHSPlasH
	Utilities
	tinyexpr
	MD5
	Ext_NeighborhoodSearch
	Ext_GenericParameters
	Ext_Discregrid
)

target_include_directories(SPlisHSPlasH PUBLIC
	${EIGEN3_INCLUDE_DIR}
	${NEIGHBORHOOD_SEARCH_INCLUDE_DIR}
	${GENERICPARAMETERS_INCLUDE_DIR}
	${DISCREGRID_INCLUDE_DIR}
)

target_link_libraries(SPlisHSPlasH
	PUBLIC Utilities tinyexpr
	INTERFACE
	${NEIGHBORHOOD_SEARCH_LIBRARIES}
	${DISCREGRID_LIBRARIES}
)

source_group("Header Files\\WCSPH" FILES ${WCSPH_HEADER_FILES})
source_group("Source Files\\WCSPH" FILES ${WCSPH_SOURCE_FILES})
source_group("Header Files\\PCISPH" FILES ${PCISPH_HEADER_FILES})
source_group("Source Files\\PCISPH" FILES ${PCISPH_SOURCE_FILES})
source_group("Header Files\\PBF" FILES ${PBF_HEADER_FILES})
source_group("Source Files\\PBF" FILES ${PBF_SOURCE_FILES})
source_group("Header Files\\IISPH" FILES ${IISPH_HEADER_FILES})
source_group("Source Files\\IISPH" FILES ${IISPH_SOURCE_FILES})
source_group("Header Files\\DFSPH" FILES ${DFSPH_HEADER_FILES})
source_group("Source Files\\DFSPH" FILES ${DFSPH_SOURCE_FILES})
source_group("Header Files\\PF" FILES ${PF_HEADER_FILES})
source_group("Source Files\\PF" FILES ${PF_SOURCE_FILES})
source_group("Header Files\\ICSPH" FILES ${ICSPH_HEADER_FILES})
source_group("Source Files\\ICSPH" FILES ${ICSPH_SOURCE_FILES})
source_group("Header Files\\Viscosity" FILES ${VISCOSITY_HEADER_FILES})
source_group("Source Files\\Viscosity" FILES ${VISCOSITY_SOURCE_FILES})
source_group("Header Files\\SurfaceTension" FILES ${SURFACETENSION_HEADER_FILES})
source_group("Source Files\\SurfaceTension" FILES ${SURFACETENSION_SOURCE_FILES})
source_group("Header Files\\Vorticity" FILES ${VORTICITY_HEADER_FILES})
source_group("Source Files\\Vorticity" FILES ${VORTICITY_SOURCE_FILES})
source_group("Header Files\\Drag" FILES ${DRAG_HEADER_FILES})
source_group("Source Files\\Drag" FILES ${DRAG_SOURCE_FILES})
source_group("Header Files\\Elasticity" FILES ${ELASTICITY_HEADER_FILES})
source_group("Source Files\\Elasticity" FILES ${ELASTICITY_SOURCE_FILES})
source_group("Header Files\\Utils" FILES ${UTILS_HEADER_FILES})
source_group("Source Files\\Utils" FILES ${UTILS_SOURCE_FILES})

