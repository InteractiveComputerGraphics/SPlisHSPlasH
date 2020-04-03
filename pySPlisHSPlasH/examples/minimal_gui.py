import pysplash as sph


def main():
    base = sph.Exec.SimulatorBase()
    base.init()
    gui = sph.GUI.Simulator_GUI_TweakBar(base)
    base.setGui(gui)
    base.run()


if __name__ == "__main__":
    main()
