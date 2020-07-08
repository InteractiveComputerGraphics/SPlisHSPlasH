import pysplishsplash as sph
import argparse
from argparse import RawTextHelpFormatter

def main():

	parser = argparse.ArgumentParser(description='SPlisHSPlasH', formatter_class=RawTextHelpFormatter)
	parser.add_argument('scene', type=str, help='scene file')
	parser.add_argument('--no-cache', action='store_true', help=
	'Disable caching of boundary samples/maps.')
	parser.add_argument('--no-initial-pause', action='store_true', help=
	'Disable initial pause when starting the simulation.')
	parser.add_argument('--no-gui', action='store_true', help=
	'Disable GUI.')
	parser.add_argument('--stopAt', type=float, default=-1.0, help='Sets or overwrites the stopAt parameter of the scene.')
	parser.add_argument('--state-file', type=str, default='', help='State file (state_<time>.bin) that should be loaded.')
	parser.add_argument('--output-dir', type=str, default='', help='Output directory for log file and partio files.')
	parser.add_argument('--param', type=str, default='', help='Sets or overwrites a parameter of the scene.\n\n' 
					  '- Setting a fluid parameter:\n\t<fluid-id>:<parameter-name>:<value>\n'
					  '- Example: --param Fluid:viscosity:0.01\n\n'
					  '- Setting a configuration parameter:\n\t<parameter-name>:<value>\n'
					  '- Example: --param cflMethod:1\n')	

	args = parser.parse_args()
	sceneFile = args.scene

	base = sph.Exec.SimulatorBase()
	
	base.init(sceneFile=args.scene, useGui=not args.no_gui, initialPause=not args.no_initial_pause, useCache=not args.no_cache, stopAt=args.stopAt, stateFile=args.state_file, outputDir=args.output_dir, param=args.param)
	
	gui = sph.GUI.Simulator_GUI_imgui(base)
	base.setGui(gui)
	base.run()


if __name__ == "__main__":
	main()
