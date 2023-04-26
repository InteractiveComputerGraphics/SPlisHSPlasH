<img src="https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/logo.jpg" width="250">
<br>

<p align=center><img src="https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/workflows/build-linux/badge.svg">&nbsp;&nbsp; <img src="https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/workflows/build-windows/badge.svg">&nbsp;&nbsp; <a href='https://splishsplash.readthedocs.io/en/latest/?badge=latest'><img src='https://readthedocs.org/projects/splishsplash/badge/?version=latest' alt='Documentation Status' /></a></p>
<p align=center>
 <img src="https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/teaser.gif">
</p>

SPlisHSPlasH is an open-source library for the physically-based simulation of fluids. The simulation in this library is based on the Smoothed Particle Hydrodynamics (SPH) method which is a popular meshless Lagrangian approach to simulate complex fluid effects. The SPH formalism allows an efficient computation of a certain quantity of a fluid particle by considering only a finite set of neighboring particles. One of the most important research topics in the field of SPH methods is the simulation of incompressible fluids. SPlisHSPlasH implements current state-of-the-art pressure solvers (WCSPH, PCISPH, PBF, IISPH, DFSPH, PF) to simulate incompressibility. Moreover, the library provides different methods to simulate viscosity, surface tension and vorticity. 

The library uses the following external libraries: [Eigen](http://eigen.tuxfamily.org/), [json](https://github.com/nlohmann/json/), [partio](https://github.com/wdas/partio/), [zlib](https://github.com/madler/zlib), [cxxopts](https://github.com/jarro2783/cxxopts), [tinyexpr](https://github.com/codeplea/tinyexpr), [toojpeg](https://github.com/stbrumme/toojpeg), [pybind](https://github.com/pybind/pybind11), [glfw](https://www.glfw.org/), [hapPLY](https://github.com/nmwsharp/happly), [nfd](https://github.com/btzy/nativefiledialog-extended), and [imgui](https://github.com/ocornut/imgui). All external dependencies are included. 

Furthermore we use our own libraries:
- [PositionBasedDynamics](https://github.com/InteractiveComputerGraphics/PositionBasedDynamics/) to simulate dynamic rigid bodies
- [Discregrid](https://github.com/InteractiveComputerGraphics/Discregrid) to detect collisions between rigid bodies
- [CompactNSearch](https://github.com/InteractiveComputerGraphics/CompactNSearch) to perform the neighborhood search 
- [cuNSearch](https://github.com/InteractiveComputerGraphics/cuNSearch) to perform the neighborhood search on the GPU
- [GenericParameters](https://github.com/InteractiveComputerGraphics/GenericParameters) to handle generic parameters

SPlisHSPlasH can export the particle data in the partio and vtk format. If you want to import partio files in Maya or Blender, try out our plugins: 
- [Blender Sequence Loader](https://github.com/InteractiveComputerGraphics/blender-sequence-loader)
- [MayaPartioTools](https://github.com/InteractiveComputerGraphics/MayaPartioTools)

**Author**: [Jan Bender](http://www.interactive-graphics.de)

## License

The SPlisHSPlasH library code is licensed under the MIT license. See [LICENSE](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/blob/master/LICENSE) for details.

External dependencies are covered by separate licensing terms.
See the [extern](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/tree/master/extern) folder for the code and respective licensing terms of each dependency.


## Documentation

* [Documentation](https://splishsplash.readthedocs.io)
* [SPH tutorial](https://interactivecomputergraphics.github.io/SPH-Tutorial)

## Forum

On our [GitHub discussions](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/discussions) page you can ask questions, discuss about simulation topics, and share ideas.


## Build Instructions

This project is based on [CMake](https://cmake.org/). Simply generate project, Makefiles, etc. using [CMake](https://cmake.org/) and compile the project with a compiler of your choice that supports C++11. The code was tested with the following configurations:
- Windows 10 64-bit, CMake 3.18.3, Visual Studio 2019
- Debian 11.5 64-bit, CMake 3.18.4, GCC 10.2.1.

Note: Please use a 64-bit target on a 64-bit operating system. 32-bit builds on a 64-bit OS are not supported.

## Python Installation Instruction

For Windows and Linux targets there exists prebuilt python wheel files which can be installed using
```
pip install pysplishsplash
```
These are available for Python Versions. See also here: [pySPlisHSPlasH](https://pypi.org/project/pySPlisHSPlasH/).
If you do not meet these conditions please refer to the build instructions and to the python binding 
[Getting started guide](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/blob/master/doc/py_getting_started.md).

The command line simulator is available by running one of the following
```
splash
splash --help
```

## Features

SPlisHSPlasH implements:
* an open-source SPH fluid simulation (2D & 3D)
* neighborhood search on CPU or GPU
* supports vectorization using AVX
* Python binding (thanks to Stefan Jeske)
* supports embedded Python scripts
* several implicit pressure solvers (WCSPH, PCISPH, PBF, IISPH, DFSPH, PF)
* explicit and implicit viscosity methods
* current surface tension approaches
* different vorticity methods
* computation of drag forces
* support for multi-phase simulations
* simulation of deformable solids 
* rigid-fluid coupling with static and dynamic bodies
* two-way coupling with deformable solids
* XSPH velocity filter
* fluid emitters
* scripted animation fields
* a json-based scene file importer
* automatic surface sampling
* a tool for volume sampling of closed geometries
* a tool to generate spray, foam and bubble particles in a postprocessing step 
* a tool to skin a visual mesh to the moving particles of an elastic solid in a postprocessing step
* partio file export of all particle data
* VTK file export of all particle data (enables the data import in ParaView)
* rigid body export
* a Maya plugin to model and generate scene files 
* a ParaView plugin to import particle data

## Pressure Solvers

The SPlisHSPlasH library implements the following pressure solvers: 

* Weakly compressible SPH for free surface flows (WCSPH)
* Predictive-corrective incompressible SPH (PCISPH)
* Position based fluids (PBF)
* Implicit incompressible SPH (IISPH)
* Divergence-free smoothed particle hydrodynamics (DFSPH)
* Projective Fluids (PF)
* Implicit compressible SPH (ICSPH)

## Boundary Handling 

The SPlisHSPlasH library implements the following boundary handling methods:

* Nadir Akinci, Markus Ihmsen, Gizem Akinci, Barbara Solenthaler, and Matthias Teschner, "Versatile rigid-fluid coupling for incompressible SPH", ACM Transactions on Graphics 31(4), 2012
* Dan Koschier and Jan Bender, "Density Maps for Improved SPH Boundary Handling", In Proceedings of ACM SIGGRAPH / EUROGRAPHICS Symposium on Computer Animation (SCA), 2017
* Jan Bender, Tassilo Kugelstadt, Marcel Weiler, Dan Koschier, "Volume Maps: An Implicit Boundary Representation for SPH", ACM SIGGRAPH Conference on Motion, Interaction and Games, 2019

## Viscosity

The SPlisHSPlasH library implements explicit viscosity methods:

* Standard SPH formulation of viscosity

and the implicit methods of the following publications:  

* Jan Bender and Dan Koschier, "Divergence-free SPH for incompressible and viscous fluids", IEEE Transactions on Visualization and Computer Graphics, 2017
* Andreas Peer, Markus Ihmsen, Jens Cornelis, and Matthias Teschner, "An Implicit Viscosity Formulation for SPH Fluids", ACM Transactions on Graphics, 34(4), 2015
* Andreas Peer and Matthias Teschner. Prescribed Velocity Gradients for Highly Viscous SPH Fluids with Vorticity Diffusion. IEEE Transactions on Visualization and Computer Graphics, 2016
* An improved version of: Tetsuya Takahashi, Yoshinori Dobashi, Issei Fujishiro, Tomoyuki Nishita, and Ming C. Lin. Implicit Formulation for SPH-based Viscous Fluids. Computer Graphics Forum, 34, 2015.
* Marcel Weiler, Dan Koschier, Magnus Brand and Jan Bender. A Physically Consistent Implicit Viscosity Solver for SPH Fluids. Computer Graphics Forum (Eurographics), 37(2), 2018

## Surface Tension

The SPlisHSPlasH library implements the surface tension methods of the following publications: 

* Markus Becker and Matthias Teschner. Weakly compressible SPH for free surface flows. In Proceedings of ACM SIGGRAPH/Eurographics Symposium on Computer Animation, 2007. Eurographics Association.
* Nadir Akinci, Gizem Akinci, and Matthias Teschner. Versatile surface tension and adhesion for SPH fluids. ACM Trans. Graph., 32(6):182:1–182:8, 2013. 
* Xiaowei He, Huamin Wang, Fengjun Zhang, Hongan Wang, Guoping Wang, and Kun Zhou, "Robust simulation of sparsely sampled thin features in SPH-based free surface flows", ACM Transactions on Graphics, 34(1), 2014.
* F. Zorilla, M. Ritter, J. Sappl, W. Rauch, M. Harders, "Accelerating  Surface Tension Calculation in SPH via Particle Classification and Monte Carlo Integration", Computers 9, 23, 2020.

## Vorticity

The SPlisHSPlasH library implements the vorticity methods of the following publications: 

* Jan Bender, Dan Koschier, Tassilo Kugelstadt and Marcel Weiler. A Micropolar Material Model for Turbulent SPH Fluids. In Proceedings of ACM SIGGRAPH / EUROGRAPHICS Symposium on Computer Animation, 2017
* Miles Macklin and Matthias Müller. Position based fluids. ACM Trans. Graph., 32(4):104:1–104:12, July 2013.

## Drag Forces

The SPlisHSPlasH library implements the drag force computation of the following publications: 

* Christoph Gissler, Stefan Band, Andreas Peer, Markus Ihmsen and Matthias Teschner. Approximate Air-Fluid Interactions for SPH. In Proceedings of Virtual Reality Interactions and Physical Simulations, 2017
* Miles Macklin, Matthias Müller, Nuttapong Chentanez and Tae-Yong Kim. Unified Particle Physics for Real-Time Applications. ACM Trans. Graph., 33(4), 2014

## Elastic Forces

* M. Becker, M. Ihmsen, and M. Teschner. Corotated SPH for deformable solids. Proceedings of Eurographics Conference on Natural Phenomena, 2009
* A. Peer, C. Gissler, S. Band, and M. Teschner. An Implicit SPH Formulation for Incompressible Linearly Elastic Solids. Computer Graphics Forum, 2017
* Tassilo Kugelstadt, Jan Bender, José Antonio Fernández-Fernández, Stefan Rhys  Jeske, Fabian Löschner, and Andreas Longva. Fast Corotated Elastic SPH  Solids with Implicit Zero-Energy Mode Control. Proceedings of the ACM  on Computer Graphics and Interactive Techniques, 2021


## Multi-Phase Fluid Simulation

The SPlisHSPlasH library implements the following publication to realize multi-phase simulations: 

* B. Solenthaler and R. Pajarola. Density Contrast SPH Interfaces. In Proceedings of ACM SIGGRAPH/Eurographics Symposium on Computer Animation, 2008.

## Volume Sampling

The SPlisHSPlasH library implements the volume sampling techniques of following publications: 

 * M. Jiang, Y. Zhou, R. Wang, R. Southern, J. J. Zhang. Blue noise sampling using an SPH-based method. ACM Transactions on Graphics, 2015
 * Tassilo Kugelstadt, Jan Bender, José Antonio Fernández-Fernández, Stefan Rhys  Jeske, Fabian Löschner, and Andreas Longva. Fast Corotated Elastic SPH  Solids with Implicit Zero-Energy Mode Control. Proceedings of the ACM  on Computer Graphics and Interactive Techniques, 2021


## Screenshots

|![](https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/SPlisHSPlasH2.jpg)|![](https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/SPlisHSPlasH1.jpg)|
|--|--|
|![](https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/SPlisHSPlasH3.jpg)|![](https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/SPlisHSPlasH4.jpg)|

## Videos

The following videos were generated using the SPlisHSPlasH library:

*A Micropolar Material Model for Turbulent SPH Fluids* | *Density Maps for Improved SPH Boundary Handling*
:---:|:---:
[![Video](https://img.youtube.com/vi/fsvDbzEui3w/0.jpg)](https://www.youtube.com/watch?v=fsvDbzEui3w) | [![Video](https://img.youtube.com/vi/P82qmTAahg0/0.jpg)](https://www.youtube.com/watch?v=P82qmTAahg0)
*Divergence-Free Smoothed Particle Hydrodynamics* | *Divergence-Free SPH for Incompressible and Viscous Fluids*
[![Video](https://img.youtube.com/vi/POnmzzhc5E0/0.jpg)](https://www.youtube.com/watch?v=POnmzzhc5E0) | [![Video](https://img.youtube.com/vi/tl4mx0TtaAc/0.jpg)](https://www.youtube.com/watch?v=tl4mx0TtaAc)
*A Physically Consistent Implicit Viscosity Solver for SPH Fluids* | *Turbulent Micropolar SPH Fluids with Foam*
[![Video](https://img.youtube.com/vi/D_nEhix1G-w/0.jpg)](https://www.youtube.com/watch?v=D_nEhix1G-w) | [![Video](https://img.youtube.com/vi/elZieJNBYqk/0.jpg)](https://www.youtube.com/watch?v=elZieJNBYqk)
*Volume Maps: An Implicit Boundary Representation for SPH* | *Implicit Frictional Boundary Handling for SPH*
[![Video](https://img.youtube.com/vi/AV_pl1bMIb8/0.jpg)](https://www.youtube.com/watch?v=AV_pl1bMIb8) | [![Video](https://img.youtube.com/vi/1u5N0eedzic/0.jpg)](https://www.youtube.com/watch?v=1u5N0eedzic) 
*Fast Corotated Elastic SPH Solids with Implicit Zero-Energy Mode Control* | 
[![Video](https://img.youtube.com/vi/8NkyiftmDN0/0.jpg)](https://www.youtube.com/watch?v=8NkyiftmDN0) | 


## References

* Nadir Akinci, Gizem Akinci, and Matthias Teschner. Versatile surface tension and adhesion for SPH fluids. ACM Trans. Graph., 32(6):182:1–182:8, 2013. 
* Nadir Akinci, Markus Ihmsen, Gizem Akinci, Barbara Solenthaler, and Matthias Teschner, "Versatile rigid-fluid coupling for incompressible SPH", ACM Transactions on Graphics 31(4), 2012
* Markus Becker and Matthias Teschner. Weakly compressible SPH for free surface flows. In Proceedings of ACM SIGGRAPH/Eurographics Symposium on Computer Animation, 2007. Eurographics Association.
* M. Becker, M. Ihmsen, and M. Teschner. Corotated SPH for deformable solids. Proceedings of Eurographics Conference on Natural Phenomena, 2009
* Jan Bender and Dan Koschier. Divergence-free smoothed particle hydrodynamics. In Proceedings of ACM SIGGRAPH / Eurographics Symposium on Computer Animation, 2015. ACM.
* Jan Bender and Dan Koschier. Divergence-free SPH for incompressible and viscous fluids. IEEE Transactions on Visualization and Computer Graphics, 2017.
* Jan Bender, Dan Koschier, Tassilo Kugelstadt and Marcel Weiler. A Micropolar Material Model for Turbulent SPH Fluids. In Proceedings of ACM SIGGRAPH / EUROGRAPHICS Symposium on Computer Animation, 2017
* Jan Bender, Dan Koschier, Tassilo Kugelstadt and Marcel Weiler. Turbulent Micropolar SPH Fluids with Foam. IEEE Transactions on Visualization and Computer Graphics 25(6), 2019
* Jan Bender, Tassilo Kugelstadt, Marcel Weiler, Dan Koschier, "Volume Maps: An Implicit Boundary Representation for SPH", ACM SIGGRAPH Conference on Motion, Interaction and Games, 2019
* Jan Bender, Tassilo Kugelstadt, Marcel Weiler, Dan Koschier, "Implicit  Frictional Boundary Handling for SPH", IEEE Transactions on  Visualization and Computer Graphics, 2020
* Jan Bender, Matthias Müller, Miguel A. Otaduy, Matthias Teschner, and Miles Macklin. A survey on position-based simulation methods in computer graphics. Computer Graphics Forum, 33(6):228–251, 2014.
* Jan Bender, Matthias Müller, and Miles Macklin. Position-based simulation methods in computer graphics. In EUROGRAPHICS 2015 Tutorials. Eurographics Association, 2015.
* Christoph Gissler, Stefan Band, Andreas Peer, Markus Ihmsen and Matthias Teschner. Approximate Air-Fluid Interactions for SPH. In Proceedings of Virtual Reality Interactions and Physical Simulations, 2017
* C. Gissler, A. Henne, S. Band, A. Peer and M. Teschner. An Implicit Compressible SPH Solver for Snow Simulation, ACM Transactions on Graphics 39(4), 2020. 
* Xiaowei He, Huamin Wang, Fengjun Zhang, Hongan Wang, Guoping Wang, and Kun Zhou. Robust simulation of sparsely sampled thin features in SPH-based free surface flows. ACM Trans. Graph., 34(1):7:1–7:9, December 2014. 
* Markus Ihmsen, Nadir Akinci, Gizem Akinci, Matthias Teschner. Unified spray, foam and air bubbles for particle-based fluids. The Visual Computer 28(6), 2012
* Markus Ihmsen, Jens Cornelis, Barbara Solenthaler, Christopher Horvath, and Matthias Teschner. Implicit incompressible SPH. IEEE Transactions on Visualization and Computer Graphics, 20(3):426–435, March 2014.
* Markus Ihmsen, Jens Orthmann, Barbara Solenthaler, Andreas Kolb, and Matthias Teschner. SPH Fluids in Computer Graphics. In Eurographics 2014 - State of the Art Reports. The Eurographics Association, 2014. 
* M. Jiang, Y. Zhou, R. Wang, R. Southern, J. J. Zhang. Blue noise sampling using an SPH-based method. ACM Transactions on Graphics, 2015
* Dan Koschier and Jan Bender, "Density Maps for Improved SPH Boundary Handling", In Proceedings of ACM SIGGRAPH / EUROGRAPHICS Symposium on Computer Animation (SCA), 2017
* Tassilo Kugelstadt, Jan Bender, José Antonio Fernández-Fernández, Stefan Rhys  Jeske, Fabian Löschner, and Andreas Longva. Fast Corotated Elastic SPH  Solids with Implicit Zero-Energy Mode Control. Proceedings of the ACM  on Computer Graphics and Interactive Techniques, 2021
* Miles Macklin and Matthias Müller. Position based fluids. ACM Trans. Graph., 32(4):104:1–104:12, July 2013.
* Miles Macklin, Matthias Müller, Nuttapong Chentanez and Tae-Yong Kim. Unified Particle Physics for Real-Time Applications. ACM Trans. Graph., 33(4), 2014
* J. J. Monaghan. Smoothed Particle Hydrodynamics. Annual Review of Astronomy and Astrophysics, 1992, 30, 543-574. 
* A. Peer, C. Gissler, S. Band, and M. Teschner. An Implicit SPH Formulation for Incompressible Linearly Elastic Solids. Computer Graphics Forum, 2017
* Andreas Peer, Markus Ihmsen, Jens Cornelis, and Matthias Teschner. An Implicit Viscosity Formulation for SPH Fluids. ACM Trans. Graph., 34(4), 2015.
* Andreas Peer and Matthias Teschner. Prescribed Velocity Gradients for Highly Viscous SPH Fluids with Vorticity Diffusion. IEEE Transactions on Visualization and Computer Graphics, 2016.
* B. Solenthaler and R. Pajarola. Density Contrast SPH Interfaces. In Proceedings of ACM SIGGRAPH/Eurographics Symposium on Computer Animation, 2008.
* B. Solenthaler and R. Pajarola. Predictive-corrective incompressible SPH. ACM Trans. Graph., 28(3):40:1–40:6, July 2009. 
* Tetsuya Takahashi, Yoshinori Dobashi, Issei Fujishiro, Tomoyuki Nishita, and Ming C. Lin. Implicit Formulation for SPH-based Viscous Fluids. Computer Graphics Forum, 34, 2015.
* Marcel Weiler, Dan Koschier and Jan Bender. Projective Fluids. Proceedings of the 9th International Conference on Motion in Games, ACM, 2016, 79-84
* Marcel Weiler, Dan Koschier, Magnus Brand and Jan Bender. A Physically Consistent Implicit Viscosity Solver for SPH Fluids. Computer Graphics Forum (Eurographics), 37(2), 2018
* F. Zorilla, M. Ritter, J. Sappl, W. Rauch, M. Harders. Accelerating  Surface Tension Calculation in SPH via Particle Classification and Monte Carlo Integration. Computers 9, 23, 2020.


## Other research projects using SPlisHSPlasH

* Diogo Schaffer, Andre Antonitsch, Amyr Neto, Soraia Musse. Towards Animating Virtual Humans in Flooded Environments. Motion, Interaction and Games, 2020
https://dl.acm.org/doi/10.1145/3424636.3426900
* Byungsoo Kim, Vinicius C. Azevedo, Markus Gross, Barbara Solenthaler. Lagrangian neural style transfer for fluids. ACM Transactions on Graphics 39, 4, 2020
https://dl.acm.org/doi/abs/10.1145/3386569.3392473
* Fernando Zorilla, Marcel Ritter, Johannes Sappl, Wolfgang Rauch, Matthias Harders. Accelerating Surface Tension Calculation in SPH via Particle Classification and Monte Carlo Integration. Computer Graphics and Visual Computing (CGVC), 2019
https://diglib.eg.org/handle/10.2312/cgvc20191260
*  H. R. Abbasia and R. Lubbad. A numerical model for the simulation of oil–ice interaction. Physics of Fluids 33, 2021
https://aip.scitation.org/doi/10.1063/5.0065587
* Uljad Berdica, Yuewei Fu, Yuchen Liu, Emmanouil Angelidis, Chen Feng. Mobile 3D Printing Robot Simulation with Viscoelastic Fluids. IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 2021
https://ieeexplore.ieee.org/document/9636114
* Arnaud Schoentgen, Pierre Poulin, Emmanuelle Darles, Philippe Meseure. Particle-based liquid control using animation templates. ACM SIGGRAPH/Eurographics Symposium on Computer Animation 2020
https://dl.acm.org/doi/10.1111/cgf.14103
* B. Ummenhofer, L. Prantl, N. Thuerey, V. Koltun. Lagrangian Fluid Simulation with Continuous Convolutions. ICLR 2020
https://ge.in.tum.de/publications/2020-ummenhofer-iclr/
*  Stefan Reinhardt, Tim Krake, Bernhard Eberhardt, Daniel Weiskopf. Consistent Shepard Interpolation for SPH-Based Fluid Animation. ACM Transactions on Graphics 38, 6, 2019
https://dl.acm.org/doi/10.1145/3355089.3356503
* Zhongyao Yang, Maolin Wu, Shiguang Liu. Helmholtz decomposition-based SPH. Virtual Reality & Intelligent Hardware 3, 2, 2021
https://www.sciencedirect.com/science/article/pii/S2096579621000176
* Min Li, Hongshu Li, Weiliang Meng, Jian Zhu, Gary Zhang. An efficient non-iterative smoothed particle hydrodynamics fluid simulation method with variable smoothing length. Visual Computing for Industry, Biomedicine, and Art 6, 1, 2023
https://vciba.springeropen.com/articles/10.1186/s42492-022-00128-x
* Yalmar Ponce Atencio, Manuel J. Ibarra, Juan José Oré Cerrón, Roberto Quispe Quispe, Richard Flores Condori, Julio Huanca Marín, Mary Huaman Carrión. Particle-Based Physics for Interactive Applications. Lecture Notes in Networks and Systems book series (LNNS,volume 216), 2021
https://link.springer.com/chapter/10.1007/978-981-16-1781-2_38
* Zijie Li, Tianqin Li, Amir Barati Farimani. TPU-GAN: Learning temporal coherence from dynamic point cloud sequences. ICLR 2022
https://openreview.net/forum?id=FEBFJ98FKx
* Muzaffer Akbay, Nicholas Nobles, Victor Zordan, Tamar Shinar. An extended partitioned method for conservative solid-fluid coupling. ACM Transactions on Graphics 37, 4, 2018
https://dl.acm.org/doi/10.1145/3197517.3201345
