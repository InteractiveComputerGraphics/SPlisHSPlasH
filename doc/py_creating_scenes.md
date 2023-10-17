# Creating Scenes

## Loading the empty scene

Right now the easiest way to create a custom scene without specifying a `Scene.json` file, is to load the predefined empty scene.

```python
import pysplishsplash as sph
import pysplishsplash.Utilities.SceneLoaderStructs as Scenes

base = sph.Exec.SimulatorBase()
base.init(sceneFile=Scenes.Empty)
```

This scene will set the default simulation method to be `DFSPH` and some other default values, which can all be changed later on.

## Recreating the double dam break scenario

In order to recreate the double dam break scenario, we need to add a bounding box as well as two fluid cubes.
The bounding box can be added as follows

```python
scene = base.getScene()
scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 3.0, 0.], scale=[4., 6., 4.], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25]))
```

The two fluid blocks can at the end be added using 

```python
scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', box=Scenes.Box([-1.5, 0.0, -1.5], [-0.5, 2.0, -0.5]), mode=0, initialVelocity=[0.0, 0.0, 0.0]))
scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', box=Scenes.Box([0.5, 0.0, 0.5], [1.5, 2.0, 1.5]), mode=0, initialVelocity=[0.0, 0.0, 0.0]))
```

This will recreate a somewhat larger scene than the default double dam break.


## Putting it all together

The following shows a script detailing how to build and run a custom double dam break.
Follow the instruction from before to activate/ deactivate the GUI.

```python
import pysplishsplash as sph
import pysplishsplash.Utilities.SceneLoaderStructs as Scenes

def main():
    # Set up the simulator
    base = sph.Exec.SimulatorBase()
    base.init(useGui=True,  sceneFile=sph.Extras.Scenes.Empty)

    # Create a simulator
    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)

    # Get the scene and add objects
    scene = base.getScene()
    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 3.0, 0.], scale=[4., 6., 4.], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25]))
    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', box=Scenes.Box([-1.5, 0.0, -1.5], [-0.5, 2.0, -0.5]), mode=0, initialVelocity=[0.0, 0.0, 0.0]))
    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', box=Scenes.Box([0.5, 0.0, 0.5], [1.5, 2.0, 1.5]), mode=0, initialVelocity=[0.0, 0.0, 0.0]))

    # Run the GUI
    base.run()

if __name__ == "__main__":
    main()
```

## Changing parameters of the simulation

The command ```base.run()``` at the end of the example above is a shortcut for the commands:
```python
   base.initSimulation()       # init simulation, fluid models, time step, ...
   base.runSimulation()        # start the simulation
   base.cleanup()              # cleanup everything after the simulation
```

If you want to change the parameters of the fluid model or the simulation, you first have to initialize the scene. Therefore, exchange ```base.run()``` by the three commands above. Then, you can change the parameters after the call of ```base.initSimulation()```.

Here is an example:

```python 
    base.initSimulation()

    # change rest density of first fluid model
	fluid = sim.getFluidModel(0)
	fluid.setValueFloat(fluid.DENSITY0, 800.0)

    # change viscosity coefficient
    visco = fluid.getViscosityBase()
	visco.setValueFloat(visco.VISCOSITY_COEFFICIENT, 0.015)

    base.runSimulation()
    base.cleanup()
```



## Loading a scene from  file

Loading a scene from a file is as simple as simply specifying a custom scene file in the init function.
This must be an **absolute path**!

```python
custom_scene = os.path.abspath("scene.json")
base.init(sceneFile=custom_scene)
```

If you want to use a gui to locate the scene file you may want to use tkinter

```python
import tkinter as tk
from tkinter import filedialog

tk.Tk().withdraw() # Dont show main window
custom_scene = filedialog.askopenfilename()
base.init(sceneFile=custom_scene)
```