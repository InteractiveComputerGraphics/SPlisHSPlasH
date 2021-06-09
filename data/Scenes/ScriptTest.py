import splishsplash as sph
import numpy as np

# List of functions that can be executed via the SPlisHSPlasH GUI (buttons will be added).
function_list = ['getExePath', 'getOutputPath', 'getSceneFile']

counter = 0
simulationBase = None

# This function is called automatically when the script is loaded in SPlisHSPlasH.
# The parameter base contains the current SimulationBase object.
def init(base):
    global counter
    global simulationBase
    
    print("init")
    simulationBase = base
    counter = 1
    
# This function is called automatically in each simulation step.     
def step():
    global counter
    sim = sph.Simulation.getCurrent()
    fluid = sim.getFluidModel(0)
    tm = sph.TimeManager.getCurrent()
    print(fluid.getPosition(0))
    
    print(counter)
    counter += 1
    print("---")
    
# This function is called automatically in each simulation reset.     
def reset():
    print("reset")
    
def getExePath():
    global simulationBase
    print("Exe path:")
    print(simulationBase.getExePath())
    print("-")
        
def getOutputPath():
    global simulationBase
    print("Output path:")
    print(simulationBase.getOutputPath())
    print("-")
    
def getSceneFile():
    print("Scene file:")
    print(sph.Exec.SceneConfiguration.getCurrent().getSceneFile())
    print("-")