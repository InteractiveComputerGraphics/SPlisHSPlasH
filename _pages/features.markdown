---
layout: splash
classes: wide
title:  "Features"
permalink: /features/
header:
  overlay_color: "#5e616c"
  overlay_image: /assets/SPlisHSPlasH.jpg
seo_desc: "Features - SPlisHSPlasH - SPH simulation of fluids and solids"
excerpt: >
  <br />
  <br />
  <br />
seo_title: My Descriptive and Keyword-Rich SEO Title  
#toc: true
#toc_label: "Table of Contents"
#toc_icon: "cog"
---
## General

* open-source SPH fluid simulation (2D & 3D)
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
* Jeske, Stefan Rhys, Lukas Westhofen, Fabian Löschner, José Antonio Fernández-Fernández, and Jan Bender. "Implicit Surface Tension for SPH Fluid Simulation." ACM Transactions on Graphics, 2023. https://doi.org/10.1145/3631936.

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
