---
permalink: /
layout: splash
hidden: true
header:
  overlay_color: "#5e616c"
  overlay_image: /assets/SPlisHSPlasH.jpg
  actions:
    - label: "<i class='fab fa-fw fa-github'></i> GitHub"
      url: "https://github.com/InteractiveComputerGraphics/SPlisHSPlasH"
excerpt: >
  An open-source library for the <br/>simulation of fluids and solids<br />
  <small><a href="https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/releases/tag/2.12.0">Latest release v2.12.0</a></small>
---
SPlisHSPlasH is an open-source library for the physically-based simulation of fluids and solids. The simulation in this library is based on the Smoothed Particle Hydrodynamics (SPH) method which is a popular meshless Lagrangian approach to simulate complex fluid effects. The SPH formalism allows an efficient computation of a certain quantity of a fluid particle by considering only a finite set of neighboring particles. One of the most important research topics in the field of SPH methods is the simulation of incompressible fluids. SPlisHSPlasH implements current state-of-the-art pressure solvers (WCSPH, PCISPH, PBF, IISPH, DFSPH, PF) to simulate incompressibility. Moreover, the library provides different methods to simulate viscosity, surface tension, vorticity, elasticity, and drag forces. 

<p align=center>
 <img src="https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/teaser.gif">
</p>

SPlisHSPlasH can export the particle data in the partio and vtk format. If you want to import partio files in Maya or Blender, try out our plugins: 
- [Blender Sequence Loader](https://github.com/InteractiveComputerGraphics/blender-sequence-loader)
- [MayaPartioTools](https://github.com/InteractiveComputerGraphics/MayaPartioTools)

A surface reconstruction of the particle data can be performed by: 
- [splashsurf](http://splashsurf.physics-simulation.org/)

SPlisHSPlasH uses the following external libraries: [Eigen](http://eigen.tuxfamily.org/), [json](https://github.com/nlohmann/json/), [partio](https://github.com/wdas/partio/), [zlib](https://github.com/madler/zlib), [cxxopts](https://github.com/jarro2783/cxxopts), [tinyexpr](https://github.com/codeplea/tinyexpr), [toojpeg](https://github.com/stbrumme/toojpeg), [pybind](https://github.com/pybind/pybind11), [glfw](https://www.glfw.org/), [hapPLY](https://github.com/nmwsharp/happly), [nfd](https://github.com/btzy/nativefiledialog-extended), and [imgui](https://github.com/ocornut/imgui). All external dependencies are included. 

Furthermore we use our own libraries:
- [PositionBasedDynamics](https://github.com/InteractiveComputerGraphics/PositionBasedDynamics/) to simulate dynamic rigid bodies
- [Discregrid](https://github.com/InteractiveComputerGraphics/Discregrid) to detect collisions between rigid bodies
- [CompactNSearch](https://github.com/InteractiveComputerGraphics/CompactNSearch) to perform the neighborhood search 
- [cuNSearch](https://github.com/InteractiveComputerGraphics/cuNSearch) to perform the neighborhood search on the GPU
- [GenericParameters](https://github.com/InteractiveComputerGraphics/GenericParameters) to handle generic parameters

**Author**: [Jan Bender](https://animation.rwth-aachen.de/person/1/)

## License

The SPlisHSPlasH library code is licensed under the MIT license. See [LICENSE](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/blob/master/LICENSE) for details.

External dependencies are covered by separate licensing terms.
See the [extern](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/tree/master/extern) folder for the code and respective licensing terms of each dependency.
