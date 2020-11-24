# This example adds a scalar and a vector particle field 
# that can be used for visualization and export
 
import pysplishsplash as sph
import numpy as np

scalar_field = np.array([])
vector_field = np.array([])
base = 0

def time_step_callback():
    global base
    global scalar_field
    global vector_field    
    sim = sph.Simulation.getCurrent()
    fluid = sim.getFluidModel(0)
    tm = sph.TimeManager.getCurrent()
    
    vel = np.array(fluid.getFieldBuffer("velocity"), copy=False)
    for i in range(0,fluid.numActiveParticles()):
        scalar_field[i] = sim.numberOfNeighbors(0, 0, i)         # output number of neighbors
        vector_field[i,:] = vel[i,:]                             # output a vector field
       
    # recompute min/max of scalar field and set corresponding values for rendering
    base.determineMinMaxOfScalarField()
    
def main():
    global base
    global scalar_field
    global vector_field
    
    base = sph.Exec.SimulatorBase()
    base.init()
    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)

    base.setTimeStepCB(time_step_callback)

    # init the simulation
    base.initSimulation()

    sim = sph.Simulation.getCurrent()
    
    # in this example we have just one fluid phase
    fluid = sim.getFluidModel(0)
    
    # resize and init field values
    # important: use the correct type: single or double 
    # depends on how you compiled SPlisHSPlasH, default is single    
    scalar_field = np.zeros(fluid.numParticles()).astype(np.single)
    vector_field = np.zeros((fluid.numParticles(), 3)).astype(np.single)
    
    # add field description
    fieldDesc = sph.FieldDescription("my_scalar_field", sph.FieldType.Scalar, sph.makeVoidPointerFct(scalar_field), False)
    fluid.addField(fieldDesc);
    fieldDesc = sph.FieldDescription("my_vector_field", sph.FieldType.Vector3, sph.makeVoidPointerFct(vector_field), False)
    fluid.addField(fieldDesc);
    
    base.setColorField(0, "my_scalar_field")
    
    gui.initSimulationParameterGUI()
    
    base.runSimulation()
    
    fluid.removeFieldByName("my_scalar_field")
    fluid.removeFieldByName("my_vector_field")
    
    base.cleanup()    

if __name__ == "__main__":
    main()
