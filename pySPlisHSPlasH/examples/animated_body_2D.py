import pysplishsplash as sph
import pysplishsplash.Utilities.SceneLoaderStructs as Scenes
import numpy as np
import math
from scipy.spatial.transform import Rotation as R


def time_step_callback():  
    sim = sph.Simulation.getCurrent()
    boundary = sim.getBoundaryModel(1)
    animatedBody = boundary.getRigidBodyObject()

    tm = sph.TimeManager.getCurrent()
    t = tm.getTime()
    #animatedBody.setAngularVelocity([0, 0, 10.0/180.0*math.pi])
    animatedBody.setVelocity([ 0.35*math.sin(0.75*t), -0.35*math.cos(0.75*t), 0])
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
    scene.particleRadius = 0.025
    scene.sim2D = True

    # change camera position
    base.setVec3ValueReal(base.CAMERA_POSITION, [0,2,8])
    base.setVec3ValueReal(base.CAMERA_LOOKAT, [0,2,0])
    base.setValueInt(base.RENDER_WALLS, 1)

    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 3.0, 0.], scale=[3., 6., 3.], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25]))

    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/Dragon_50k.obj", translation=[-0.45, 1, 0.], scale=[1.5,1.5,1.5], color=[0.5, 0.5, 0.5, 1.0], isAnimated=True, isWall=False, mapInvert=False, mapResolution=[25, 25, 25]))

    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin = [-1.5, 0.0, -0.1], boxMax = [1.5 , 1.0, 0.1], mode=0, initialVelocity=[0.0, 0.0, 0.0]))

    scene.emitters.append(Scenes.EmitterData(id='Fluid', x=[1,3,0], width=5, height=5, axis = [0,1,0], angle = math.pi, emitStartTime=0.0, emitEndTime=3.0, velocity=3))

    scene.materials.append(Scenes.MaterialData(id='Fluid', colorMapType=1, maxVal=8.0))

    # init the simulation
    base.initSimulation()

    sim = sph.Simulation.getCurrent()
    sim.setValueInt(sim.BOUNDARY_HANDLING_METHOD, 0)

    base.runSimulation()
    base.cleanup()

if __name__ == "__main__":
    main()

