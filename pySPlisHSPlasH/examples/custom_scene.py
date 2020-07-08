import pysplishsplash as sph
import pysplishsplash.Utilities.SceneLoaderStructs as Scenes


def main():
    # Set up the simulator
    base = sph.Exec.SimulatorBase()
    base.init(useGui=True,  sceneFile=sph.Extras.Scenes.Empty)

    # Create an imgui simulator
    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)

    # Get the scene and add objects
    scene = base.getScene()
    scene.boundaryModels.append(Scenes.BoundaryData(meshFile="../models/UnitBox.obj", translation=[0., 3.0, 0.], scale=[4., 6., 4.], color=[0.1, 0.4, 0.5, 1.0], isWall=True, mapInvert=True, mapResolution=[25, 25, 25]))
    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', box=Scenes.Box([-1.5, 0.0, -1.5], [-0.5, 2.0, -0.5]), mode=0, initialVelocity=[0.0, 0.0, 0.0]))
    scene.fluidBlocks.append(Scenes.FluidBlock(id='Fluid', box=Scenes.Box([0.5, 0.0, 0.5], [1.5, 2.0, 1.5]), mode=0, initialVelocity=[0.0, 0.0, 0.0]))

    # Run the GUI
    base.run()


if __name__ == "__main__":
    main()
