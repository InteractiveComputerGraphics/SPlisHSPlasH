# CMake Options

This page should give you a short overview over the CMake options of SPlisHSPlasH.

## USE_DOUBLE_PRECISION

If this flag is enabled, then all computations with floating point values are performed using double precision (double). Otherwise single precision (float) is used.

## USE_AVX

SPlishSPlasH supports the usage of AVX (Advanced Vector Extensions) which is an extension of modern CPUs to perform a single instruction on multiple data. The extension allows to perform eight floating point operations in parallel. Enabling AVX significantly improves the performance of the simulator. Currently, the following methods have AVS support: 

* DFSPH 
* the micropolar vorticity model
* the standard viscosity model
* the viscosity model of Weiler et al.

## USE_OpenMP

Enable the OpenMP parallelization which lets the simulation run in parallel on all available cores of the CPU. 

## USE_GPU_NEIGHBORHOOD_SEARCH

As default SPlisHSPlasH uses [CompactNSearch ](https://github.com/InteractiveComputerGraphics/CompactNSearch) as neighborhood search which performs all operations on the CPU. However, with this flag you can switch to [cuNSearch ](https://github.com/InteractiveComputerGraphics/cuNSearch) which is our GPU neighborhood search. In case you want to use the GPU method, you have to install Cuda.

## USE_IMGUI

We just reimplemented the GUI using [imgui](https://github.com/ocornut/imgui) instead of [AntTweakBar](http://anttweakbar.sourceforge.net). If you want to try out the new GUI, enable this flag. 
