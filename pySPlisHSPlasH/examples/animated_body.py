import pysplishsplash as sph
import pysplishsplash.Utilities.SceneLoaderStructs as Scenes
import numpy as np
import math
from scipy.spatial.transform import Rotation as R

def time_step_callback():  
    sim = sph.Simulation.getCurrent()
    boundary = sim.getBoundaryModel(1)
    animatedBody = boundary.getRigidBodyObject()
   
    animatedBody.setAngularVelocity([0, 60.0/180.0*math.pi, 0])
    animatedBody.animate()

def main():
    # Set up the simulator
    base = sph.Exec.SimulatorBase()
    base.init(useGui=True,  sceneFile=sph.Extras.Scenes.Empty)

    # Create an imgui simulator
    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)
    base.setTimeStepCB(time_step_callback)

    # Get the scene and add objects
    scene = sph.Exec.SceneConfiguration.getCurrent().getScene()
    
    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 3.0, 0.], scale=[3., 6., 3.], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25]))
    r = R.from_euler('y', 30.0, degrees=True).as_matrix()
    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/Dragon_50k.obj", translation=[0., 0.05, 0.], scale=[1,1,1], color=[0.5, 0.5, 0.5, 1.0], rotation=r, isAnimated=True, isWall=False, mapInvert=False, mapResolution=[25, 25, 25]))
    
    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', box=Scenes.Box([-1.25, 0.0, -1.25], [-0.25, 1.5, -0.25]), mode=0, initialVelocity=[0.0, 0.0, 0.0]))

    # init the simulation
    base.initSimulation()

    sim = sph.Simulation.getCurrent()
    sim.setValueInt(sim.BOUNDARY_HANDLING_METHOD, 2)
    
    base.runSimulation()
    base.cleanup()

if __name__ == "__main__":
    main()

