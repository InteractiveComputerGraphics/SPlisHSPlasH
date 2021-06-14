# Getting started

This page should give you a short overview of SPlisHSPlasH.

SPlisHSPlasH currently consists of a simulators and different tools which are introduced in the following:

## SPHSimulator

This application reads a SPlisHSPlasH scene file and performs a simulation of the scene. 

The scene file format is explained [here.](file_format.md)

##### Command line options:

* -h, --help: Print help text.
* -v, --version: Print version.
* --no-cache: Disable caching of boundary samples/maps.
* --state-file: Load a simulation state of the corresponding scene.
* --output-dir: Output directory for log file and partio files.
* --no-initial-pause: Disable caching of boundary samples/maps.
* --no-gui: Disable graphical user interface. The simulation is run only in the command line without graphical output. The "stopAt" option must be set in the scene file or by the next parameter.
* --stopAt arg: Sets or overwrites the stopAt parameter of the scene.
* --param arg: Sets or overwrites a parameter of the scene.
	- Setting a fluid parameter: <fluid-id>:<parameter-name>:<value>
		- Example: --param Fluid:viscosity:0.01
	- Setting a configuration parameter: <parameter-name>:<value>
		- Example: --param cflMethod:1

##### Hotkeys

* Space: pause/contiunue simulation
* r: reset simulation
* w: wireframe rendering of meshes
* m: recompute min and max values for color-coding the color field in the rendering process
* i: print all field information of the selected particles to the console
* s: save current simulation state
* l: load simulation state (currently only Windows)
* +: perform a single time step
* ESC: exit

## Python bindings 

SPlisHSPlasH implements bindings for python using [pybind11](https://github.com/pybind/pybind11).
See the [getting started guide](./py_getting_started.md).

### Impatient installation guide

In order to install, simply clone the repository and run pip install on the repository.
It is recommended, that you set up a **virtual environment** for this, because cache files will be stored in the directory of the python installation along with models and scene files.

```shell script
git clone https://github.com/InteractiveComputerGraphics/SPlisHSPlasH.git
pip install SPlisHSPlasH/
```





