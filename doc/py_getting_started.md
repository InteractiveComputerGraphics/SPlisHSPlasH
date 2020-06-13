# pySPlisHSPlasH

## Python bindings for the SPlisHSPlasH library

## Requirements

Currently the generation of python bindings is only tested on 

- Linux Debian, gcc 8.3, Python 3.7/3.8 (Anaconda), CMake 3.13
- Windows 10, Visual Studio 15/17/19, Python 3.7/3.8 (Anaconda), CMake 3.13

Note that the compiler, the python installation as well as cmake have to be available from the command
line for the installation process to work. 
MacOS builds should work but have not been tested.

## Installation

In order to install it is advised that you create a new virtual environment so that any faults during
installation can not mess up your python installation.
This is done as follows for 

**conda**

```shell script
conda create --name venv python=3.7
conda activate venv
```

**virtualenv**

```shell script
python3 -m virtualenv venv --python=python3.7
source venv/bin/activate
```

Now you can clone the repository by

```shell script
git clone https://github.com/InteractiveComputerGraphics/SPlisHSPlasH.git
```

And finally you should be able to install SPlisHSPlasH using pip. 
**The trailing slash is important** otherwise pip will try to download the package, which is not supported yet at least.
Also note, that `pip install SPlisHSPlasH` should be called from **one directory above** the cloned source directory and **not within** the directory itself.

```shell script
pip install SPlisHSPlasH/
```

While `pip install` is useful if SPlisHSPlasH should only be installed once, for development purposes it might be more sensible to build differently.
Change into the SPlisHSPlasH directory and build a python wheel file as follows

```shell script
cd SPlisHSPlasH
python setup.py bdist_wheel
pip install -I build/dist/*.whl
```

When building a new version of SPlisHSPlasH simply run these commands again and the installation will be updated.
The compile times will be lower, because the build files from previous installations remain. 
If you are getting compile errors please try to compile the pysplishsplash target of the CMake project separately.

Now check your installation by running

```shell script
python -c "import pysplishsplash"
```

**Note**: You may have to install numpy. 
Future releases may already contain numpy as a dependency.

```shell script
pip install numpy
```

## I want to see something very very quickly

If you're very impatient, just run the following command after installing

```shell script
splash
```

You will be prompted to select a preconfigured scene file which will then be run in a User Interface.
For more options and functionality run.
The keybindings in the GUI are the same as for the regular SPlisHSPlasH version.

```shell script
splash --help
```

## Minimal working example

The following examples should work, if SPlisHSPlasH was installed correctly. 
If you want to load other scene files, be sure to place them into the SPlisHSPlasH data directory structure.

**With GUI**

```python
import pysplishsplash as sph

def main():
    base = sph.Exec.SimulatorBase()
    base.init()
    gui = sph.GUI.Simulator_GUI_TweakBar(base)
    base.setGui(gui)
    base.run()

if __name__ == "__main__":
    main()
```

**Without GUI**

```python
import pysplishsplash as sph

def main():
    base = sph.Exec.SimulatorBase()
    base.init(useGui=False)
    base.setValueFloat(base.STOP_AT, 10.0) # Important to have the dot to denote a float
    base.run()

if __name__ == "__main__":
    main()
```

**Outputting the results to a specific directory without GUI**
```python
import pysplishsplash as sph
from pysplishsplash.Extras import Scenes
import os

def main():
    base = sph.Exec.SimulatorBase()
    output_dir = os.path.abspath("where/you/want/the/data")
    base.init(useGui=False, outputDir=output_dir, sceneFile=Scenes.DoubleDamBreak)
    base.setValueFloat(base.STOP_AT, 20.0) # Important to have the dot to denote a float
    base.setValueBool(base.VTK_EXPORT, True)
    # Uncomment the next line to set the output FPS value (must be float)
    # base.setValueFloat(base.DATA_EXPORT_FPS, 10000.) 
    base.run()

if __name__ == "__main__":
    main()
```
## SPHSimulator.py

If you want to start the simulator in the same way as the C++ version, just use the SPHSimulator.py in the examples directory.

## Modifying other properties

The bindings cover most of the public interface of the SPlisHSPlasH library.
As such, it is possible to change components of the simulation dynamically.
In the following example, the second cube in the well known double dam break scenario is replaced with a slightly larger cube.

```python
import pysplishsplash
import pysplishsplash.Utilities.SceneLoaderStructs as Scene

def main():
    base = pysplishsplash.Exec.SimulatorBase()
    args = base.init()
    gui = pysplishsplash.GUI.Simulator_GUI_TweakBar(base)
    base.setGui(gui)
    scene = base.getScene()
    add_block = Scene.FluidBlock('Fluid', Scene.Box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), 0, [0.0, 0.0, 0.0])
    scene.fluidBlocks[1] = add_block # In Place construction not supported yet
    base.run()

if __name__ == "__main__":
    main()
```