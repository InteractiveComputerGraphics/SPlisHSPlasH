import splishsplash as sph
import numpy as np
import math

simulationBase = None

# This function is called automatically when the script is loaded in SPlisHSPlasH.
# The parameter base contains the current SimulationBase object.
def init(base):
    global simulationBase
    simulationBase = base
    
# This function is called automatically in each simulation step.     
def step():
    sim = sph.Simulation.getCurrent()
    boundary = sim.getBoundaryModel(1)
    animatedBody = boundary.getRigidBodyObject()
   
    tm = sph.TimeManager.getCurrent()
    t = tm.getTime()
    animatedBody.setVelocity([ 0.35*math.sin(0.75*t), -0.35*math.cos(0.75*t), 0])
    animatedBody.animate()
