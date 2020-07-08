import pysplishsplash as sph

def key_callback():
	print("Hello World")
	
def time_step_callback():
	print("step")

def main():
	base = sph.Exec.SimulatorBase()
	base.init()
	gui = sph.GUI.Simulator_GUI_imgui(base)
	base.setGui(gui)
	
	gui.addKeyFunc('k', key_callback)
	base.setTimeStepCB(time_step_callback)
	
	base.run()


if __name__ == "__main__":
	main()
