# Getting started

This page should give you a short overview of SPlisHSPlasH.

SPlisHSPlasH currently consists of a simulators and different tools which are introduced in the following:

## SPHSimulator

This application reads a SPlisHSPlasH scene file and performs a simulation of the scene. 

The scene file format is explained [here.](file_format.md)

##### Command line options:

* -h, --help: Print help text.

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
* i: print all field information of the selected particles to the console
* s: save current simulation state
* l: load simulation state (currently only Windows)
* ESC: exit

## Python bindings 

SPlisHSPlasH implements bindings for python using [pybind11](https://github.com/pybind/pybind11).
See the [getting started guide](./pysplash/getting_started.md).

### Impatient installation guide

In order to install, simply clone the repository and run pip install on the repository.
It is recommended, that you set up a **virtual environment** for this, because cache files will be stored in the directory of the python installation along with models and scene files.

```shell script
git clone https://github.com/InteractiveComputerGraphics/SPlisHSPlasH.git
pip install SPlisHSPlasH/
```



## Tools

## partio2vtk

A tool to convert partion files in vtk files. In this way the particle data which is exported from SPlisHSPlasH can be converted to the vtk format. This is useful to import the data in ParaView for visualization.

## PartioViewer

The simulators can export the particle simulation data using the partio file format. The PartioViewer can read such a file and render the particle data using OpenGL. This tool is able to handle multiphase data and rigid body data. It can create image sequences and movies (using ffmpeg).

To visualize a sequence of partio files or a single file, call (the index in the file name is used for the sequence): 
```
PartioViewer fluid_data_1.bgeo
```

This tool is also able to read a complete output directory:
```
PartioViewer output/DamBreakModel
```
In this case the tool searches for the partio files of multiple phases in the subdirectory "partio" and for rigid body data in "rigid_bodies".

Note: To generate videos you must tell PartioViewer where it can find the ffmpeg executable.

##### Command line options:

* -h, --help: Print help
* --renderSequence: Render a sequence from startFrame to endFrame as jpeg.
* --renderVideo: Render a sequence from startFrame to endFrame as video.This function requires ffmpeg which must be in the PATH or the ffmpegPath parameter must be set.
* --noOverwrite: Do not overwrite existing frames when using --renderSequence option. Existing frames are not loaded at all which accelerates the image sequence generation.
* -o, --outdir arg: Output directory for images
* --rbData arg: Rigid body data to visualize (bin file)
* --ffmpegPath arg: Path of the ffmpeg excutable.
* --width arg: Width of the image in pixels. (default: 1024)
* --height arg: Height of the image in pixels. (default: 768)
* --fps arg: Frame rate of video. (default: 25)
* -r, --radius arg: Particle radius (default: 0.025)
* -s, --startFrame arg: Start frame (only used if value is >= 0) (default: -1)
* -e, --endFrame arg: End frame (only used if value is >= 0) (default: -1)
* --colorField arg: Name of field that is used for the color. (default: velocity)
* --colorMapType arg: Color map (0=None, 1=Jet, 2=Plasma) (default: 1)
* --renderMinValue arg: Min value of field. (default: 0.0)
* --renderMaxValue arg: Max value of field. (default: 10.0)
* --camPos arg: Camera position (e.g. --camPos "0 1 5") (default: 0 3 10)
* --camLookat arg: Camera lookat (e.g. --camLookat "0 0 0") (default: 0 0 0)

##### Hotkeys

* Space: pause/contiunue simulation
* r: reset simulation
* w: wireframe rendering of meshes
* i: print all field information of the selected particles to the console
* s: save current frame as jpg image
* v: generate video 
* j: generate image sequence
* +: step to next frame
* -: step to previous frame
* ESC: exit


## SurfaceSampling

A popular boundary handling method which is also implemented in SPlisHSPlasH uses a particle sampling of the surfaces of all boundary objects. This command line tool can generate such a surface sampling. Note that the same surface sampling is also integrated in the simulators and the samplings are generated automatically if they are required. However, if you want to generate a surface sampling manually, then you can use this tool. 

## VolumeSampling

The simulators can load particle data from partio files. This particle data then defines the initial configuration of the particles in the simulation. The VolumeSampling tool allows you to sample a volumetric object with particle data. This means you can load an OBJ file with a closed surface geometry and sample the interior with particles. 