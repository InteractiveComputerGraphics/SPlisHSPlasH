pybind11_add_module(pysplishsplash
  main.cpp
  AnimationFieldModule.cpp
  BoundaryModelModule.cpp
  ParameterObjectModule.cpp
  UtilitiesModule.cpp
  SPHKernelsModule.cpp
  EmitterModule.cpp
  FluidModelModule.cpp
  NonPressureForceBaseModule.cpp
  SimulationModule.cpp
  StaticRigidBodyModule.cpp
  TimeModule.cpp
  TriangleMeshModule.cpp
  DFSPHModule.cpp
  DragModule.cpp
  ElasticityModule.cpp
  IISPHModule.cpp
  PBFModule.cpp
  PCISPHModule.cpp
  PFModule.cpp
  SurfaceTensionModule.cpp
  ViscosityModule.cpp
  VorticityModule.cpp
  WCSPHModule.cpp
  GUIModule.cpp
  ExtrasModule.cpp
)

# All the same postfix, otherwise the module name ist not the same as the exported name and python gets confused
set_target_properties(pysplishsplash PROPERTIES
        DEBUG_POSTFIX ""
        RELWITHDEBINFO_POSTFIX ""
        MINSIZEREL_POSTFIX "")

target_link_libraries(pysplishsplash PRIVATE SPlisHSPlasH SimulatorBase)
