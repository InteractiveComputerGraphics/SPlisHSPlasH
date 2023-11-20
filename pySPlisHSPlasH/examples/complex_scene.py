import pysplishsplash as sph
import pysplishsplash.Utilities.SceneLoaderStructs as Scenes
import numpy as np
import math
from scipy.spatial.transform import Rotation as R


def main():
	# Set up the simulator
	base = sph.Exec.SimulatorBase()
	base.init(useGui=True, sceneFile=sph.Extras.Scenes.Empty)
	
	sim = sph.Simulation.getCurrent()

	# Create an imgui simulator
	gui = sph.GUI.Simulator_GUI_imgui(base)
	base.setGui(gui)

	# SimulationBase
	base.setValueBool(base.PAUSE, False)
	base.setValueFloat(base.PAUSE_AT, -1.0)
	base.setValueFloat(base.STOP_AT, 10.0)
	base.setValueUInt(base.NUM_STEPS_PER_RENDER, 4)
	base.setValueInt(base.RENDER_WALLS, 4)
	base.setValueFloat(base.DATA_EXPORT_FPS, 25.0)
	base.setValueBool(base.STATE_EXPORT, False)
	base.setValueFloat(base.STATE_EXPORT_FPS, 25.0)
	base.setValueString(base.PARTICLE_EXPORT_ATTRIBUTES, "velocity;density")
	base.activateExporter("Partio Exporter", False)
	base.activateExporter("VTK Exporter", False)
	base.activateExporter("Rigid Body Exporter", False)
	base.activateExporter("Rigid Body OBJ Exporter", False)	
	base.activateExporter("Rigid Body VTK Exporter", False)	

	# Get the scene and add objects
	scene = sph.Exec.SceneConfiguration.getCurrent().getScene()
	
	# change camera position
	base.setVec3ValueReal(base.CAMERA_POSITION, [2,3,6])
	base.setVec3ValueReal(base.CAMERA_LOOKAT, [0,0,0])
	
	# change time step size
	tm = sph.TimeManager.getCurrent()
	tm.setValueFloat(tm.TIME_STEP_SIZE, 0.002)

	# change parameters of scene
	scene.particleRadius = 0.025
	scene.sim2D = False
	
	# add objects
	scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 3.0, 0.], scale=[4., 6., 4.], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25], isDynamic=False))

	# first fluid
	scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', boxMin=[-1.5, 0.0, -1.5], boxMax = [-0.5, 1.0, -0.5], mode=0, initialVelocity=[0.0, 0.0, 0.0]	))
	scene.fluidModels.append(Scenes.FluidData(id='Fluid', samplesFile="../models/sphere.obj", mode=0, scale=[0.3, 0.3, 0.3], translation=[0., 0.5, 0.], initialVelocity=[0.0, 0.0, 0.0]))
	
	scene.materials.append(Scenes.MaterialData(id='Fluid', colorMapType=2))
	
	# second fluid (emitter)	
	scene.emitters.append(Scenes.EmitterData(id='Fluid2', x=[1,1,0], width=5, height=5, axis=[0,1,0], angle=math.pi, emitStartTime=0.0, emitEndTime=1.0, velocity=3))
	
	scene.materials.append(Scenes.MaterialData(id='Fluid2', colorMapType=1, minVal=0.0, maxVal=10.0, maxEmitterParticles=5000))
	
	# add animation field
	scene.animatedFields.append(Scenes.AnimationFieldData(particleFieldName='velocity', startTime=0, endTime=5, expressionX='cos(2*t)'))
	scene.animatedFields[0].scale=[2,2,2]
	scene.animatedFields[0].translation=[-1,0,0]
	
	# init the simulation
    # This creates the fluid models, timeStep, ... 
    # Call base.initSimulation() before sim.getTimeStep() or sim.getFluidModel()
	base.initSimulation()
		
	# Simulation
	sim.setValueBool(sim.ENABLE_Z_SORT, False)	
	sim.setValueInt(sim.CFL_METHOD, 2)
	sim.setValueFloat(sim.CFL_FACTOR, 0.5)
	sim.setValueFloat(sim.CFL_MIN_TIMESTEPSIZE, 0.0002)
	sim.setValueFloat(sim.CFL_MAX_TIMESTEPSIZE, 0.005)
	sim.setValueInt(sim.KERNEL_METHOD, 4)
	sim.setValueInt(sim.GRAD_KERNEL_METHOD, 4)
	sim.setValueInt(sim.SIMULATION_METHOD, 4)
	sim.setValueInt(sim.BOUNDARY_HANDLING_METHOD, 2)

	# DFSPH
	timeStep = sim.getTimeStep()
	timeStep.setValueUInt(timeStep.MIN_ITERATIONS, 2)
	timeStep.setValueUInt(timeStep.MAX_ITERATIONS, 50)
	timeStep.setValueFloat(timeStep.MAX_ERROR, 0.05)
	timeStep.setValueUInt(timeStep.MAX_ITERATIONS_V	, 50)
	timeStep.setValueFloat(timeStep.MAX_ERROR_V, 0.1)
	timeStep.setValueBool(timeStep.USE_DIVERGENCE_SOLVER, True)
	
	# first FluidModel 
	fluid = sim.getFluidModel(0)
	fluid.setValueFloat(fluid.DENSITY0, 800.0)
	fluid.setValueInt(fluid.DRAG_METHOD, 0)
	fluid.setValueInt(fluid.SURFACE_TENSION_METHOD, 0)
	fluid.setValueInt(fluid.VORTICITY_METHOD, 1)
	fluid.setValueInt(fluid.ELASTICITY_METHOD, 0)
	fluid.setViscosityMethod("Standard")
	
	# Viscosity
	visco = fluid.getViscosityBase()
	visco.setValueFloat(visco.VISCOSITY_COEFFICIENT, 0.015)
	
	# Vorticity
	vorticity = fluid.getVorticityBase()
	vorticity.setValueFloat(vorticity.VORTICITY_COEFFICIENT, 0.05)
	
	# second FluidModel 
	fluid2 = sim.getFluidModel(1)
	fluid2.setValueFloat(fluid2.DENSITY0, 1000.0)
	fluid2.setValueInt(fluid2.VISCOSITY_METHOD, 1)
	
	# Viscosity
	visco2 = fluid2.getViscosityBase()
	visco2.setValueFloat(visco2.VISCOSITY_COEFFICIENT, 0.01)
	
	base.runSimulation()
	base.cleanup()


if __name__ == "__main__":
    main()
