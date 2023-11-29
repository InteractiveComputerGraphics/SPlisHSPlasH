---
layout: single
classes: wide
title:  "Getting started"
permalink: /getting_started/
layout: splash
header:
  overlay_color: "#5e616c"
  overlay_image: /assets/SPlisHSPlasH.jpg
seo_desc: "Getting started - SPlisHSPlasH - SPH simulation of fluids and solids"  
excerpt: >
  <br />
  <br />
  <br />
#toc: true
#toc_label: "Table of Contents"
#toc_icon: "cog"
---
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
These are available for Python versions 3.6-3.10. See also here: [pySPlisHSPlasH](https://pypi.org/project/pySPlisHSPlasH/).
If you do not meet these conditions please refer to the build instructions and to the python binding 
[Getting started guide](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH/blob/master/doc/py_getting_started.md).

The command line simulator is available by running one of the following
```
splash
splash --help
```