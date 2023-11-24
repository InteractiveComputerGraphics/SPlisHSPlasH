# Replicability

The SPlisHSPlasH library implements the SPH methods developed by our and other research groups (build instructions can be found [here](build_from_source.md)). This allows to reproduce the research results of the corresponding publications. Inspired by the [Graphics Replicability Stamp Initiative](http://www.replicabilitystamp.org/) we started to add scenes to the repository to reproduce some of the results in our papers. 

The easiest way to run these scenes is by using the precompiled python packages `pip install pysplishsplash` and running `splash` from the command line, as described in [quick-start](../README.md#python-installation-instruction) or the [python documentation](py_getting_started.md)

**Jan Bender, Tassilo Kugelstadt, Marcel Weiler, Dan Koschier, "Implicit Frictional Boundary Handling for SPH", IEEE Transactions on Visualization and Computer Graphics, 2020**

* Figure 7.a) can be replicated by loading the scene: data/Scenes/GridModel_Akinci2012.json
* Figure 7.b) can be replicated by loading the scene: data/Scenes/GridModel_Bender2019.json

**Jeske, Stefan Rhys, Lukas Westhofen, Fabian Löschner, José Antonio Fernández-Fernández, and Jan Bender. “Implicit Surface Tension for SPH Fluid Simulation.” ACM Transactions on Graphics, 2023**

* Figure 8 (left) can be replicated by loading the scene: `data/Scenes/SurfaceTension_CoveredSphere_JWL+23.json`
* Figure 9 can be replicated by loading the scene: `data/Scenes/SurfaceTension_WaterBell_JWL+23.json`
* Figure 10 can be replicated by loading the scene: `data/Scenes/SurfaceTension_FluidChain_JWL+23.json`
