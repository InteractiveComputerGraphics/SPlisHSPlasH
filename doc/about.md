![](https://raw.githubusercontent.com/InteractiveComputerGraphics/SPlisHSPlasH/master/doc/images/logo_small.jpg)

# About SPlisHSPlasH

SPlisHSPlasH is an open-source library for the physically-based simulation of fluids. The simulation in this library is based on the Smoothed Particle Hydrodynamics (SPH) method which is a popular meshless Lagrangian approach to simulate complex fluid effects. The SPH formalism allows an efficient computation of a certain quantity of a fluid particle by considering only a finite set of neighboring particles. One of the most important research topics in the field of SPH methods is the simulation of incompressible fluids. SPlisHSPlasH implements current state-of-the-art pressure solvers (WCSPH, PCISPH, PBF, IISPH, DFSPH, PF) to simulate incompressibility. Moreover, the library provides different methods to simulate viscosity, surface tension and vorticity. 

## Main features

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
* fluid emitters
* scripted animation fields
* a json-based scene file importer
* automatic surface sampling
* a tool for volume sampling of closed geometries
* partio file export of all particle data
* VTK file export of all particle data (enables the data import in ParaView)
* rigid body export
* a Maya plugin to model and generate scene files 
* a ParaView plugin to import particle data

## License

The MIT License (MIT)

Copyright (c) 2016-present, SPlisHSPlasH contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
