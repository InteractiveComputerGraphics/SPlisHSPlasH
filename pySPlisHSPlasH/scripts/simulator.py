import pysplishsplash as sph
import sys
import os
import tkinter as tk
from tkinter import filedialog


def main():
    base = sph.Exec.SimulatorBase()

    is_windows = sys.platform.startswith("win32") or sys.platform.startswith("cygwin")
    no_gui_flags = ["-h", "--help", "--no-gui", "-v", "--version"]
    if not any(flag in sys.argv for flag in no_gui_flags) and not is_windows:
        tk.Tk().withdraw()
        scene_dir = os.path.join(os.path.dirname(sys.executable), "data/Scenes/")
        scene_file = filedialog.askopenfilename(initialdir=scene_dir)
        base.init(sys.argv + [scene_file], "[Python] SPlisHSPlasH")
    else:
        base.init(sys.argv, "[Python] SPlisHSPlasH")

    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)

    base.run()


if __name__ == "__main__":
    main()
