# Note that for this example you have to install matplotlib by:
# pip install matplotlib
 
import pysplishsplash as sph
import matplotlib.pyplot as plt
import numpy as np

counter = 0
velocities = []
times = []

def time_step_callback():
	global counter
	sim = sph.Simulation.getCurrent()
	tm = sph.TimeManager.getCurrent()
	fluid = sim.getFluidModel(0)
	
	# plot the norm of the velocity of particle 0
	v = np.array(fluid.getVelocity(0))
	vn = np.linalg.norm(v)
	velocities.append(vn)
	times.append(tm.getTime())
	
	if counter % 20 == 0:
		plt.plot(times, velocities, 'b')
		plt.draw()
		plt.pause(0.001)
	counter += 1

def main():
	plt.xlabel('time (s)')
	plt.ylabel('velocity (m/s)')
	plt.title('Example plot')
	plt.grid(True)
	plt.ion()
	plt.plot([0.0], [0.0], 'b')
	plt.draw()
	plt.pause(0.001)
	plt.show()

	base = sph.Exec.SimulatorBase()
	base.init()
	gui = sph.GUI.Simulator_GUI_imgui(base)
	base.setGui(gui)
	
	base.setTimeStepCB(time_step_callback)
	
	# init the simulation
	base.initSimulation()
	
	sim = sph.Simulation.getCurrent()
	sim.setValueBool(sim.ENABLE_Z_SORT, False)	
	
	base.runSimulation()
	base.cleanup()

if __name__ == "__main__":
	main()
