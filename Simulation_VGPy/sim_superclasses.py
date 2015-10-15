"""
Vincente Pericoli
Oct 7, 2015

Superclasses useful for abaqus simulations
"""
#
# imports
#

import numpy
import sys
import myPaths
sys.path.append(myPaths.OdbTools())
from odbFieldVariableClasses import *
from inpPartMeshClasses import *
from calcVGI import *

#
# main class
#

class superSpecimen(object):
    """ a superclass for test specimen simulations
    
    Attributes:
        abqPath      = string location of ODB file. Must be of the form:
                       'C:\\temp\\folder\\odbfile.odb'
        material     = string of material type: 'AP50' or 'AP70HP'
        setName      = string of the name of the set of interest
                       (e.g. set to obtain VGI)
    
    Attributes set by self.calcNodalMonoVGI():
        VGI          = array of entire VGI history
                       rows are associated with abaqus frames
                       columns are associated with the nodes in the set
    
    Attributes set by self.calcElementalMonoVGI():
        VGI          = same as above, but values are associated with elements
    
    Attributes set by self.fetchVolume():
        elemVol      = numpy array of element volumes in self.setName
        
    Attributes set by self.fetchMeshInfo():
        nodesCoords  = numpy array of the nodal coordinates
        elemConnect  = numpy array of the elemental connectivity
    """
    
    #
    # Attributes (object initialization)
    #
    def __init__(self):
        """ return object with desired attributes """
        
        self.abqPath  = None
        self.partName = None
        
        # must be uppercase
        self.setName  = None
        self.material = None
        
        # set from Methods:
        #calcNodalMonoVGI, calcElementalMonoVGI
        self.VGI           = None
        self.nodeLabels    = None
        self.elementLabels = None
        #fetchMesh
        self.nodesCoords   = None
        self.elemConnect   = None
        #fetchVolume
        self.elemVol       = None
        return
    
    #
    # Set Dependent Attributes
    #
    @property
    def odbName(self):
        """ returns odb file name (with file extensions) """
        return self.abqPath.split('\\')[-1]
        
    @property
    def name(self):
        """ returns name of simulation (no file extension) """
        return os.path.splitext(self.odbName)[0]

    #
    # Methods
    #
    def calcNodalMonoVGI(self):
        """ Obtains the monotonic VGI of (nodal) self.setName """
        
        # obtain the PEEQ history
        PEEQ = IntPtVariable(self.abqPath, 'PEEQ', self.setName)
        PEEQ.fetchNodalAverage()
        
        # obtain the mises history
        mises = IntPtVariable(self.abqPath, 'MISES', self.setName)
        mises.fetchNodalAverage()
        
        # obtain the pressure history
        pressure = IntPtVariable(self.abqPath, 'PRESS', self.setName)
        pressure.fetchNodalAverage()
        
        # obtain the VGI history of the simulation
        VGI = calcMonotonicVGI(mises.resultData, pressure.resultData, PEEQ.resultData)        
        
        # save VGI and labels, then return
        self.VGI = VGI
        self.nodeLabels = PEEQ.nodeLabels
        return
    
    def calcElementalMonoVGI(self):
        """ Obtains the monotonic VGI of (elemental) self.setName """
        
        # obtain the PEEQ history
        PEEQ = IntPtVariable(self.abqPath, 'PEEQ', self.setName)
        PEEQ.fetchElementAverage()
        
        # obtain the mises history
        mises = IntPtVariable(self.abqPath, 'MISES', self.setName)
        mises.fetchElementAverage()
        
        # obtain the pressure history
        pressure = IntPtVariable(self.abqPath, 'PRESS', self.setName)
        pressure.fetchElementAverage()
        
        # obtain the VGI history of the simulation
        VGI = calcMonotonicVGI(mises.resultData, pressure.resultData, PEEQ.resultData)        
        
        # save VGI and labels, then return
        self.VGI = VGI
        self.elementLabels = PEEQ.elementLabels
        return

    def fetchMeshInfo(self, partName=None):
        """ obtain the nodal coordinates and elemental connectivity """
        
        #check input
        if partName is None:
            partName = self.partName
        
        #inpPath is the name of the input file (will be same as odb)
        inpPath  = self.abqPath.replace(".odb",".inp")
        
        #generate PartMesh object and fetch the mesh
        mesh = PartMesh(inpPath, partName)
        mesh.fetchMesh()

        #save to self, return
        self.elemConnect = mesh.elemConnect
        self.nodesCoords = mesh.nodesCoords
        return
        
    def fetchVolume(self):
        """ obtain the initial volume for the elements in the self.setName """
        
        vol = ElementVariable(self.odbName, 'EVOL', self.setName)
        vol.fetchInitialElementVolume()
        self.elemVol = vol.resultData
        return
