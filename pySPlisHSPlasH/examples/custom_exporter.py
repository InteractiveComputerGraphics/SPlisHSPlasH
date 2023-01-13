import pysplishsplash as sph
import numpy as np
import meshio
import sys
import os

class CustomExporter(sph.ExporterBase):
    def __init__(self, base):
        sph.ExporterBase.__init__(self, base)
        self.exportPath = ""
        self.firstFrame = True
        
    def init(self, outputPath): 
        self.exportPath = os.path.join(outputPath, "rigid_bodies")
        
    def step(self, frame): 
        if (not self.getActive()):
            return;
        
        fileName = "rb_data_"
        exportFileName = os.path.join(self.exportPath, fileName)
        self.writeRigidBodies(exportFileName, frame)
        
    def reset(self): 
        self.firstFrame = True
        return
        
    def setActive(self, active): 
        sph.ExporterBase.setActive(self,active)
        if (self.getActive()):
            os.makedirs(self.exportPath, exist_ok = True)
            
    def writeRigidBodies(self, fileName, frame):
        sim = sph.Simulation.getCurrent()
        nBoundaryModels = sim.numberOfBoundaryModels()

        # check if we have a static model
        isStatic = True
        for i in range(0, nBoundaryModels):
            bm = sim.getBoundaryModel(i)
            if (bm.getRigidBodyObject().isDynamic()):
                isStatic = False
                break
                
           
        if (self.firstFrame or not isStatic):
            for i in range(0, nBoundaryModels):
                bm = sim.getBoundaryModel(i)
                rbo = bm.getRigidBodyObject()
                vertices = np.array(rbo.getVertexBuffer(), copy=False)
                faces = np.array(rbo.getFaceBuffer(), copy=False)
                tris = faces.reshape((-1,3))
                cells = [("triangle", tris)]
                
                meshio.write_points_cells(fileName + str(i) + "_" + str(frame) + ".obj", vertices, cells)
        self.firstFrame = False

    
def main():
    base = sph.Exec.SimulatorBase()
    
    exporter = CustomExporter(base)
    base.addRigidBodyExporter("enableCustomExporter", "Custom Exporter", "Enable/disable custom export.", exporter)

    base.init(sys.argv, "[Python] SPlisHSPlasH")

    gui = sph.GUI.Simulator_GUI_imgui(base)
    base.setGui(gui)
    
    base.activateExporter("Custom Exporter", True)	
    
    print ("Active rigid body exporters: ")
    exporters = base.getRigidBodyExporters()
    for ex in exporters:
        if ex.exporter.getActive():
            print (ex.name)
    print("")
        
    # init the simulation
    base.initSimulation()

    base.runSimulation()
    base.cleanup()

if __name__ == "__main__":
	main()
