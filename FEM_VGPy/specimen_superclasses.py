"""
Vincente Pericoli
Oct 7, 2015

Superclasses useful for abaqus simulations of structural test specimens
"""
#
# imports
#

import numpy
import sys
import myPaths
sys.path.append(myPaths.OdbTools())
from odbFieldVariableClasses import *
from odbInstanceMeshClasses import *
from calcVGI import *

#
# main class
#

class superSpecimen(object):
    """ a superclass for test specimen simulations
    
    Attributes:
        odbPath      = string location of ODB file. Must be of the form:
                       'C:\\temp\\folder\\odbfile.odb'
        material     = string of material type: 'AP50' or 'AP70HP'
        instanceName = string name of the instance in the ABAQUS assembly
        setName      = string of the name of the set of interest
                       (e.g. set to obtain VGI)
    
    Attributes set by self.calcNodalExtrapMonoVGI():
        VGI          = rank-3 array of nodal (extrapolated) VGI history 
                       for each element in the defined element set setName.
                       access is: VGI[frame,node,element]
                       see:       self.elemLabelSet, self.nodeLabelSet
                       VGI[:,:,e] corresponds to element self.elemLabelSet[e]
                       
    Attributes set by self.calcIntPtMonoVGI():
        VGI          = rank-3 array of integration point VGI history
                       for each element in the defined element set setName.
                       access is: VGI[frame,ip,element]
                       see:       self.elemLabelSet, self.intPtLabelSet
                       VGI[:,:,e] corresponds to element self.elemLabelSet[e]
                       
    Attributes set by self.calcAllMonoVGI():
        VGI          = dictionary with two keys--
                       VGI['ELEM_NODAL'] is the VGI from self.calcNodalExtrapMonoVGI
                       VGI['ELEM_IP'] is the VGI from self.calcIntPtMonoVGI
    
    Attributes set by self.calcNodalAvgMonoVGI():
        VGI          = array of entire VGI history
                       rows are associated with ABAQUS frames
                       columns are associated with the nodes in the set
    
    Attributes set by self.calcElemAvgMonoVGI():
        VGI          = same as above, but average values are associated
                       with elements
    
    Attributes set by self.fetchVolume():
        elemVol      = numpy array of element volumes in self.setName
        
    Attributes set by self.fetchMeshInfo():
        nodesCoords  = numpy array of the nodal coordinates
        elemConnect  = numpy array of the elemental connectivity
    """
    
    #
    # Attributes (object initialization)
    #
    def __init__(self, odbPath, material, instanceName,
                       setName, failureLoad=None, loadSetName=None):
        """ return object with desired attributes """
        
        self.odbPath  = odbPath
        
        # must be uppercase
        self.setName      = setName.upper()
        self.material     = material.upper()
        self.instanceName = instanceName.upper()
        if loadSetName is not None:
            self.loadSetName  = loadSetName.upper()
        
        # set from Methods:
        #calc VGI's
        self.VGI           = None
        self.failureIndex  = None
        self.nodeLabelSet  = None
        self.elemLabelSet  = None
        self.intPtLabelSet = None
        #fetchMesh
        self.nodesCoords   = None
        self.elemConnect   = None
        self.elemType      = None
        #fetchVolume
        self.elemVol       = None
        
        # initialize failureLoad related attributes.
        # these should be set/handled by the subclasses
        self.failureLoad = failureLoad
        self.loadHist    = None
        return
    
    #
    # Set Dependent Attributes
    #
    @property
    def odbName(self):
        """ returns odb file name (with file extensions) """
        return self.odbPath.split('\\')[-1]
        
    @property
    def name(self):
        """ returns name of simulation (no file extension) """
        return os.path.splitext(self.odbName)[0]

    #
    # Methods
    #
    def calcNodalExtrapMonoVGI(self):
        """ 
        Obtians an extrapolated monotonic VGI for nodes of elements 
        in (elemental) self.setName
        """
        
        # obtain the PEEQ history
        PEEQ = IntPtVariable(self.odbPath, 'PEEQ', self.setName)
        PEEQ.fetchNodalExtrap()
        
        # obtain the mises history
        mises = IntPtVariable(self.odbPath, 'MISES', self.setName)
        mises.fetchNodalExtrap()
        
        # obtain the pressure history
        pressure = IntPtVariable(self.odbPath, 'PRESS', self.setName)
        pressure.fetchNodalExtrap()
        
        # obtain the VGI history of the simulation
        VGI = calcMonotonicVGI( mises.resultData, pressure.resultData, PEEQ.resultData )
    
        # save VGI and labels, then return
        self.VGI = VGI
        self.elemLabelSet = mises.elementLabels
        self.nodeLabelSet = mises.nodeLabels
        return
        
    def calcIntPtMonoVGI(self):
        """ 
        Obtains the monotonic VGI for integration points of elements
        in (elemental) self.setName
        """
        
        # obtain the PEEQ history
        PEEQ = IntPtVariable(self.odbPath, 'PEEQ', self.setName)
        PEEQ.fetchIntPtData()
        
        # obtain the mises history
        mises = IntPtVariable(self.odbPath, 'MISES', self.setName)
        mises.fetchIntPtData()
        
        # obtain the pressure history
        pressure = IntPtVariable(self.odbPath, 'PRESS', self.setName)
        pressure.fetchIntPtData()
        
        # obtain the VGI history of the simulation
        VGI = calcMonotonicVGI( mises.resultData, pressure.resultData, PEEQ.resultData )
    
        # save VGI and labels, then return
        self.VGI = VGI
        self.elemLabelSet  = mises.elementLabels
        self.intPtLabelSet = mises.intPtLabels
        return

    def calcAllMonoVGI(self):
        """ 
        Obtains the monotonic VGI at integration points AND nodes
        (extrapolated) of elements in (elemental) self.setName
        
        i.e. it calculates the monotonic VGI at all possible 
        data locations (no averaging scheme)
        
        saves this to self.VGI as a dictionary
        integration point data has dictionary key 'ELEM_IP'
        nodal data has dictionary key 'ELEM_NODAL'
        """
        
        # obtain the nodal (extrapolated) VGI
        self.calcNodalExtrapMonoVGI()
        nodal_VGI = self.VGI
        
        # obtain the integration point VGI
        self.calcIntPtMonoVGI()
        ip_VGI = self.VGI
        
        # save as dict to VGI
        self.VGI = {'ELEM_IP':ip_VGI, 'ELEM_NODAL':nodal_VGI}
        return
        
        
    def calcNodalAvgMonoVGI(self):
        """ Obtains the average monotonic VGI of (nodal) self.setName """
        
        # obtain the PEEQ history
        PEEQ = IntPtVariable(self.odbPath, 'PEEQ', self.setName)
        PEEQ.fetchNodalAverage()
        
        # obtain the mises history
        mises = IntPtVariable(self.odbPath, 'MISES', self.setName)
        mises.fetchNodalAverage()
        
        # obtain the pressure history
        pressure = IntPtVariable(self.odbPath, 'PRESS', self.setName)
        pressure.fetchNodalAverage()
        
        # obtain the VGI history of the simulation
        VGI = calcMonotonicVGI(mises.resultData, pressure.resultData, PEEQ.resultData)        
        
        # save VGI and labels, then return
        self.VGI = VGI
        self.nodeLabelSet = PEEQ.nodeLabels
        return
    
    def calcElemAvgMonoVGI(self):
        """ Obtains the average monotonic VGI of (elemental) self.setName """
        
        # obtain the PEEQ history
        PEEQ = IntPtVariable(self.odbPath, 'PEEQ', self.setName)
        PEEQ.fetchElementAverage()
        
        # obtain the mises history
        mises = IntPtVariable(self.odbPath, 'MISES', self.setName)
        mises.fetchElementAverage()
        
        # obtain the pressure history
        pressure = IntPtVariable(self.odbPath, 'PRESS', self.setName)
        pressure.fetchElementAverage()
        
        # obtain the VGI history of the simulation
        VGI = calcMonotonicVGI(mises.resultData, pressure.resultData, PEEQ.resultData)        
        
        # save VGI and labels, then return
        self.VGI = VGI
        self.elemLabelSet = PEEQ.elementLabels
        return
        
    def fetchMeshInfo(self, instanceName=None, exactKey=False):
        """ obtain the nodal coordinates and elemental connectivity """
        # check input
        if instanceName is None:
            instanceName = self.instanceName

        # generate InstanceMesh object and fetch the mesh
        mesh = InstanceMesh(self.odbPath, instanceName, exactKey)
        mesh.fetchMesh()

        # save to self, return
        self.elemConnect = mesh.elemConnect
        self.elemType    = mesh.elemType
        self.nodesCoords = mesh.nodesCoords
        return
        
    def fetchVolume(self):
        """ obtain the initial volume for the elements in the self.setName """
        
        vol = ElementVariable(self.odbName, 'EVOL', self.setName)
        vol.fetchInitialElementVolume()
        self.elemVol = vol.resultData
        return
