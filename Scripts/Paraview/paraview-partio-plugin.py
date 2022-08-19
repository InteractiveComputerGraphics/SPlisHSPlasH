import numpy as np
import partio
import collections

from paraview.util.vtkAlgorithm import (
    VTKPythonAlgorithmBase,
    smdomain,
    smhint,
    smproperty,
    smproxy,
)
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

paraview_plugin_version = "1.0"
partio_extensions = "bgeo"

@smproxy.reader(
    name="Partio reader",
    extensions=partio_extensions,
    file_description="Partio files",
    support_reload=True,
)
class PartioReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkUnstructuredGrid"
        )
        self.fileList = []  
        self.timeSteps = None        
            
    # function is required to give Paraview the info about the time steps
    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self.timeSteps   
          
    # read the file sequence 
    @smproperty.xml(
        """
        <StringVectorProperty animateable="0"
                                clean_command="RemoveAllFileNames"
                                command="AddFileName"
                                name="FileName"
                                number_of_elements="1"
                                panel_visibility="never"
                                repeat_command="1">
            <FileListDomain name="files" />
        </StringVectorProperty>
        """ 
    )
    def AddFileName(self, name):
        self.fileList.append(name)
        self.Modified()

    def RemoveAllFileNames(self):
        self.fileList.clear()

    # file chooser in the Paraview GUI
    @smdomain.filelist()
    @smhint.filechooser(
        extensions=partio_extensions, file_description="Partio files"
    )
       
    # get current time for data
    def getCurrentTime(self, outInfo):
        executive = self.GetExecutive()
        timesteps = self.timeSteps
        if timesteps is None or len(timesteps) == 0:
            return None
        elif outInfo.Has(executive.UPDATE_TIME_STEP()) and len(timesteps) > 0:
            utime = outInfo.Get(executive.UPDATE_TIME_STEP())
            dtime = timesteps[0]
            for atime in timesteps:
                if atime > utime:
                    return dtime
                else:
                    dtime = atime
            return dtime
        else:
            assert(len(timesteps) > 0)
            return timesteps[0]
            
    # set time steps in Paraview to animate sequence
    def RequestInformation(self, request, inInfoVec, outInfoVec):
        if self.fileList is None:
            return 1
        if len(self.fileList) > 1:                      
            executive = self.GetExecutive()
            outInfo = outInfoVec.GetInformationObject(0)
            outInfo.Remove(executive.TIME_STEPS())
            outInfo.Remove(executive.TIME_RANGE())

            self.timeSteps = np.arange(0, len(self.fileList), 1, dtype=int) 
            if self.timeSteps is not None:
                for t in self.timeSteps:
                    outInfo.Append(executive.TIME_STEPS(), t)
                outInfo.Append(executive.TIME_RANGE(), self.timeSteps[0])
                outInfo.Append(executive.TIME_RANGE(), self.timeSteps[-1])
        return 1

    # transfer the data to Paraview
    def RequestData(self, request, inInfoVec, outInfoVec):
        data_time = self.getCurrentTime(outInfoVec.GetInformationObject(0))     
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfoVec))

        if (data_time != None):
            currentFile = self.fileList[data_time]
        else:
            currentFile = self.fileList[0]
        
        p = partio.read(currentFile)

        if p == None:
            return 1
            
        totalParticles = p.numParticles()
            
        for i in range(p.numAttributes()):
            attr=p.attributeInfo(i)
            if attr.name=="position":
                pos = np.array(p.data_buffer(attr), copy=False)  
                output.SetPoints(pos)
        
        # cell_conn contains tuples (num indices, vertex ids)
        # e.g. particles (1,0), (1,1), (1,2), (1,3), ...
        # e.g. triangles (3, 0, 1, 2), (3, 2, 3, 4), ...
        cell_conn = np.hstack([np.ones((totalParticles,1)), np.arange(0, totalParticles, 1, dtype=int).reshape(-1,1)]).flatten()
        # for particles use type VERTEX=1
        cell_types = np.full((totalParticles), 1, np.ubyte)
        # offset between two particles is 2 since cell_conn always contains the number of indices 
        cell_offsets = 2 * np.arange(totalParticles, dtype=int) 
        
        output.SetCells(cell_types, cell_offsets, cell_conn)
        
        # add field data
        for i in range(p.numAttributes()):
            attr=p.attributeInfo(i)
            if attr.name != "position": 
                values = np.array(p.data_buffer(attr), copy=False)
                output.PointData.append(values, attr.name)                
        return 1
   

