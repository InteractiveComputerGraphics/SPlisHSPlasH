# Getting started

This page should give you a short overview of SPlisHSPlasH.

SPlisHSPlasH currently consists of different simulators and tools which are introduced in the following:

## Simulators 

### StaticBoundarySimulator

This application reads a SPlisHSPlasH scene file and performs a simulation of the scene. It assumes that only static boundary objects are in the scenario which increases the performance. If you want to simulation dynamic boundaries, you can use "DynamicBoundarySimulator". 

The scene file format is explained [here.](file_format.md)

### DynamicBoundarySimulator

This application can also simulate SPlisHSPlasH scenes but in contrast to the StaticBoundarySimulator it can handle dynamic boundaries. The dynamic rigid bodies are simulated using our [PositionBasedDynamics library](https://github.com/InteractiveComputerGraphics/PositionBasedDynamics) which is automatically included in the build process. If a scene only contains static bodies, you should use "StaticBoundarySimulator" since it is faster. 

The scene file format is explained [here.](file_format.md)

## Tools

## PartioViewer

The simulators can export the particle simulation data using the partio file format. The PartioViewer can read such a file and render the particle data using OpenGL. 

## SurfaceSampling

A popular boundary handling method which is also implemented in SPlisHSPlasH uses a particle sampling of the surfaces of all boundary objects. This command line tool can generate such a surface sampling. Note that the same surface sampling is also integrated in the simulators and the samplings are generated automatically if they are required. However, if you want to generate a surface sampling manually, then you can use this tool. 

## VolumeSampling

The simulators can load particle data from partio files. This particle data then defines the initial configuration of the particles in the simulation. The VolumeSampling tool allows you to sample a volumetric object with particle data. This means you can load an OBJ file with a closed surface geometry and sample the interior with particles. 