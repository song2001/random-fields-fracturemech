"""
Vincente Pericoli
UC Davis


Subclasses of simulation specimens:
    * SNTT (smooth notched tensile test specimen)
    * CT   (compact tension specimen)
    * BN   (blunted notch specimen)
    * BB   (bolt-bearing specimen)
    * BH   (bolt-hole specimen)
    * RBS  (reduced beam section specimen)
"""

#
# imports
#
import numpy
import sys
import myPaths
sys.path.append(myPaths.OdbTools())
from specimen_superclasses import *
from odbFieldVariableClasses import *

#
# subclass definitions
#

class SNTT(superSpecimen):
    """ a smooth notched tensile test specimen. inherits from superSpecimen
    
    SNTT(odbPath, material, failureLoad, setName, loadSetName)
    
    Attributes:
        odbPath      = string location of ODB file. Must be of the form:
                       'C:\\temp\\folder\\odbfile.odb'
        material     = string of material type: 'AP50' or 'AP70HP'
        failureLoad  = list or tuple of the failure displacements (test results)
                       for this specimen
        setName      = string of the name of the set of interest
                       (e.g. set to obtain VGI)
        loadSetName  = string of the name of the set where displacement is applied
    """
    #
    # Attributes (object initialization)
    #
    def __init__(self, odbPath, material, failureLoad, 
                       setName='CenterNode',
                       loadSetName='DisplacementSurface'):
        """ return object with desired attributes """
        # initialize using superclass
        superSpecimen.__init__(self, odbPath=odbPath, instanceName="SNTT_", 
                                     setName=setName, material=material,
                                     failureLoad=failureLoad, loadSetName=loadSetName)
        # set an SNTT-specific radius attribute
        self.radius = float(self.name.split('_')[1][1:])/1000.0
        return
    
    #
    # Properties
    #
    @property
    def ERR_TOL(self):
        # set the failure percent error tolerance
        return 0.0005
    #
    # Methods
    #
    def fetchLoadHist(self):
        """
        fetch the abaqus displacement (load/response) history of the simulation.
        
        this is specific to SNTT specimens, such that the relevant displacement is
        taken to be in the 2-direction (AKA y-direction), and that the first node
        in the set is representative of all other nodes in the set (since they are
        all prescribed the same displacement load in ABAQUS).
        """
        
        # obtain the displacement history of the simulation.
        # assume that the relevant displacement is in the 2-direction
        # (AKA y-direction), and further assume that the first node 
        # (column) is representative of the other columns
        abqDispl = NodalVariable(self.odbPath, 'U', self.loadSetName)
        abqDispl.fetchNodalOutput()
        abqDispl = abqDispl.resultData[:,0,1]
        
        # ensure proper shape
        loadHist = numpy.zeros((max(abqDispl.shape),1),dtype=numpy.float64)
        loadHist[:,0] = abqDispl
        # save to attribute
        self.loadHist = loadHist
        return
    
    def determineFailureIndex(self):
        """
        determine which "history" (AKA frame) index corresponding to failure.
        
        these indices correspond to the nearest percent error (within tolerance)
        between the loadHist and failureLoad
        """
        
        # obtain the displacement history of the simulation, if needed
        if self.loadHist is None:
            self.fetchLoadHist()
        
        # number of frames in the abaqus history
        nframe = self.loadHist.shape[0]
        
        # determine which frames correspond to failure
        failureIndex = []
        for displ in self.failureLoad:
            # for all failure displacements, locate which frame index it corresponds to.
            # we don't want to compare floats directly, so use a percent error
            percent_err = numpy.absolute( 100.0*(self.loadHist - displ)/displ )
            
            # lowest percent error => equivalent to failure displ
            frameind = numpy.argmin(percent_err)
            
            # check to see that it is within tolerance
            if percent_err[frameind] < self.ERR_TOL:
                # good. save it.
                failureIndex.append(int(frameind))
            else:
                # provide warning that something went wrong.
                # append nothing. try to continue.
                print ("\n!! WARNING: " + self.name + ": the failure index could not " +
                       "be located for failure displacement " + str(displ) + " !!")
                print "Nearest Percent Error is: " + str(percent_err[frameind]) + "\n"
        
        # save to atttribute
        self.failureIndex = tuple(failureIndex)
        return

class CT(superSpecimen):
    """ a compact tension specimen (ASTM E1820)
    
    Attributes:
        odbPath      = string location of ODB file. Must be of the form:
                       'C:\\temp\\folder\\odbfile.odb'
        material     = string of material type: 'AP50' or 'AP70HP'
        failureLoad  = list or tuple of the failure J1 (test results)
                       for this specimen (i.e. the observed J1c's)
        setName      = string of the name of the set of interest
                       (i.e. where to obtain VGI)
        crackTipSet  = string name of the set defining the crack tip (must be a single node)
        crackName    = string name of the crack, from which to obtain the ABAQUS J1 history
        stepName     = string name of the relevant step, from which to obtain the ABAQUS J1 history
        
    failureLoad corresponds to J1 for CT specimens
    """
    #
    # Attributes (object initialization)
    #
    def __init__(self, odbPath, material, failureLoad, 
                       setName='CrackExtensionPlane', 
                       crackTipSet='CrackTip', crackName='', stepName='Pull',
                       loadSetName=None):
        """ return object with desired attributes """
        
        # initialize generic superclass variables
        superSpecimen.__init__(self, odbPath=odbPath, instanceName="CT_", 
                                     setName=setName, material=material,
                                     failureLoad=failureLoad)
        
        # initialize CT-specific variables
        self.stepName    = stepName
        self.crackName   = crackName
        self.crackTipSet = crackTipSet.upper() # must be uppercase
        # loadSetName is not meaningful for CT
        del self.loadSetName
        return

    #
    # Dependent Properties
    #
    
    @property
    def lstar(self):
        """ this is the deterministic length scale """
        if self.material == 'AP50':
            return (0.0033, 0.007, 0.017)
        elif self.material == 'AP70HP':
            return (0.0025, 0.012, 0.016)
        else:
            raise Exception('Undefined material!')
    
    #
    # Methods
    #
    def fetchDeterministicVGI(self, overwrite=True):
        """ 
        obtain the deterministic VGI.
        this means the VGI associated with the l* characteristic lengths.
        by default, this will overwrite the self.VGI and self.nodeLabels attribute
        otherwise, it will save new attributes: self.deterministicVGI
        
        this will also save a new attribute: self.lstarNodeInfo, which is [index, nodeLabel]
        """
        
        #
        # check if pre-requisites are properly met
        #
        if self.elementLabels is not None:
            #make sure element calcs have not been executed
            raise Exception("method not supported or meaningful for element output")
        elif self.nodeLabels is None:
            #if nodal calcs have not been executed, execute them...
            print "\nExecuting calcNodalAvgMonoVGI()..."
            self.calcNodalAvgMonoVGI()
            print "done!\n"
        
        #
        # determine which nodal locations are associated with the l*'s
        #
        
        # find out crack node initial coordinates
        dummy = NodalVariable(self.odbName, 'COORD', self.crackTipSet)
        dummy.fetchNodalOutput()
        crackTipCoords = dummy.resultData[0,:,:] #first frame
        del dummy
        
        # find initial coordinates of the nodes ahead of the crack tip
        dummy = NodalVariable(self.odbName, 'COORD', self.setName)
        dummy.fetchNodalOutput()
        setCoords = dummy.resultData[0,:,:] #first frame
        del dummy
        
        # find out how many l*'s there are
        nlstar = len(self.lstar)
        
        # find out which nodes are associated with the l*'s, and their index
        # (i.e., which of the nodes do we want to save data for?)
        lstarNodeInfo = numpy.zeros((2,nlstar),dtype=numpy.int_)
        
        for i,lstar in enumerate(self.lstar):
            # find out which nodes are closest to "lstar away from the crack-tip"
            # this assumes crack tip is at a smaller x-coord than the crack opening
            target = crackTipCoords[0,0] - lstar
            percent_err = numpy.absolute( 100.0*(setCoords[:,0] - target) / target )
            
            # the node we want would have the lowest percent error
            nodind = numpy.argmin(percent_err)
            
            #save the nodal index value
            lstarNodeInfo[0,i] = nodind
            #save the actual node label, for debugging...
            lstarNodeInfo[1,i] = self.nodeLabels[nodind]
        
        #
        # fetch the VGI at those node locations:
        #
        nframe = self.VGI.shape[0]
        deterministicVGI = numpy.zeros((nframe,nlstar),dtype=numpy.float64)
        
        for i,nodind in enumerate(lstarNodeInfo[0,:]):
            deterministicVGI[:,i] = self.VGI[:,nodind]
        
        #
        # save and return
        #
        if overwrite:
            del self.VGI
            self.VGI = deterministicVGI
        else:
            self.deterministicVGI = deterministicVGI
        
        self.lstarNodeInfo    = lstarNodeInfo
        return

    def fetchLoadHist(self):
        """
        fetch the abaqus J1 (load/response) history of the simulation.
        
        this is specific to CT specimens, such that the relevant load/response is
        the J1 of the defined crack.
        """
        
        # obtain the J1 history of the simulation
        cv = CrackVariable(self.odbName, self.stepName, self.crackName)
        cv.getJintegral()
        
        # number of total frames in the abaqus history
        nframeHist = cv.shape[0]
        
        # we want to prepend the J1 with a zero, since it is a history 
        # variable; history variables start at the first non-zero frame,
        # while field variables start at frame 0. This allows us to 
        # return field-variable compatible indices.
        abqJ1       = numpy.zeros((nframeHist+1,1),dtype=numpy.float64)
        abqJ1[1:,:] = cv.resultData[:,-1]

        # save to attribute
        self.loadHist = abqJ1
        return        
        # find out which frames represent observed failure
        failureIndex = []
        for J1c in self.failureLoad:
            #for each observed J1c, find the nearest abaqus J1 and it's index
            percent_err = numpy.absolute( 100.0*(abqJ1 - J1c)/J1c )
            
            # the frame we want has the lowest percent error
            index = numpy.argmin(percent_err)
            failureIndex.append( int(index) )
            
        # therefore, the VGI's at those frames are the failure VGI's... save them!
        self._save_failureVGI(tuple(failureIndex))
        
        # we also want to save the "loading" history
        self.loadHist = abqJ1
        return
        
    def determineFailureIndex(self):
        """
        determine which "history" (AKA frame) index corresponds to failure.
        
        this is specific to CT because failure is specifically 
        related to the observed J1c. this will automatically
        call self.fetchLoadHist() if it has not been already executed.
        """
        
        # obtain load history if not already set.
        if self.loadHist is None:
            self.fetchLoadHist()
        
        # rename for convenience
        abqJ1 = self.loadHist
        
        # find out which frames correspond to observed failure
        failureIndex = []
        fi_pctErr    = []
        for J1c in self.failureLoad:
            #for each observed J1c, find the nearest abaqus J1 and it's index
            percent_err = numpy.absolute( (abqJ1 - J1c)/J1c )
            
            # the frame we want has the lowest percent error
            index = numpy.argmin(percent_err)
            failureIndex.append( int(index) )
            # in the case of a CT specimen, it would be nice to save the pct error
            fi_pctErr.append( percent_err[index] )
            
        # save to attribute
        self.failureIndex       = tuple(failureIndex)
        self.failureIndexPctErr = tuple(fi_pctErr)
        return

class BN(superSpecimen):
    """ blunted notch specimen """
    #
    # Attributes (Object Initialization)
    #
    def __init__(self, odbPath, material, failureLoad,
                        setName='CenterPlane',
                        loadSetName='DisplControl'):
        """ return instance with desired attributes """
        
        # initialize using superclass with BN defaults
        superSpecimen.__init__(self, odbPath=odbPath, instanceName="BN", 
                                     setName=setName, material=material,
                                     failureLoad=failureLoad, loadSetName=loadSetName)
        return
    
    #
    # Methods
    #
    def fetchLoadHist(self):
        """
        ...
        """
        return
    
    def determineFailureIndex(self):
        """
        ...
        """
        return

class BB(SNTT):
    """ 
    bolt-bearing specimen.
    inherit from SNTT, since same determineFailureIndex() will be used.
    """
    #
    # Attributes (Object Initialization)
    #
    def __init__(self, odbPath, material, failureLoad,
                        setName='DeterministicLine',
                        loadSetName=('LVDT_top','LVDT_bottom')):
        """ return instance with desired attributes """
        
        # initialize using superclass, but with BB defaults
        superSpecimen.__init__(self, odbPath=odbPath, instanceName="BB_symm", 
                                     setName=setName, material=material,
                                     failureLoad=failureLoad)
        
        # set names must be upper-case
        self.loadSetName = tuple([ s.upper() for s in loadSetName ])        
        return
    
    #
    # Properties
    #
    @property
    def ERR_TOL(self):
        # set the failure percent error tolerance
        return 0.5
    
    #
    # Methods
    #
    def fetchLoadHist(self):
        """
        fetch the abaqus displacement (load/response) history of the simulation.
        
        this is specific to BB specimens, such that the relevant displacement is
        taken to be in the 2-direction (AKA y-direction), and the load is taken
        to be the difference between the average displacement of the top and 
        bottom LVDT sets.
        """
        
        # obtain the displacement history of the simulation.
        # assume that the relevant displacement is in the 2-direction
        # (AKA y-direction).
        lvdt0 = NodalVariable(self.odbPath, 'U', self.loadSetName[0])
        lvdt0.fetchNodalOutput()
        lvdt0.avgNodalOutput()
        lvdt0 = lvdt0.resultData[:,0,1]
        
        lvdt1 = NodalVariable(self.odbPath, 'U', self.loadSetName[1])
        lvdt1.fetchNodalOutput()
        lvdt1.avgNodalOutput()
        lvdt1 = lvdt1.resultData[:,0,1]
        
        lvdt = numpy.absolute(lvdt0 - lvdt1)
        
        # ensure proper shape
        loadHist = numpy.zeros((max(lvdt.shape),1),dtype=numpy.float64)
        loadHist[:,0] = lvdt
        # save to attribute
        self.loadHist = loadHist
        return

class BH(SNTT):
    """ 
    bolt-hole specimen.
    inherit from SNTT, since same determineFailureIndex() will be used.
    """
    #
    # Attributes (Object Initialization)
    #
    def __init__(self, odbPath, material, failureLoad,
                        setName='DeterministicLine',
                        loadSetName='LVDT'):
        """ override SNTT init. return BH instance with desired attributes """
        
        # initialize using superclass
        superSpecimen.__init__(self, odbPath=odbPath, instanceName="BH_symm", 
                                     setName=setName, material=material,
                                     failureLoad=failureLoad, loadSetName=loadSetName)
        return
        
    #
    # Properties
    #
    @property
    def ERR_TOL(self):
        # set the failure percent error tolerance
        return 0.5
    
    #
    # Methods
    #
    def fetchLoadHist(self):
        """
        override SNTT fetchLoadHist
        fetch the abaqus displacement (load/response) history of the simulation.
        
        this is specific to BB specimens, such that the relevant displacement is
        taken to be in the 2-direction (AKA y-direction), and the load is taken
        to be the 2x the displacement of the loadSet (due to symmetry)
        """
        
        # obtain the displacement history of the simulation.
        # assume that the relevant displacement is in the 2-direction
        # (AKA y-direction).
        lvdt = NodalVariable(self.odbPath, 'U', self.loadSetName)
        lvdt.fetchNodalOutput()
        lvdt.avgNodalOutput()
        lvdt = 2.0 * lvdt.resultData[:,0,1]
        
        # ensure proper shape
        loadHist = numpy.zeros((max(lvdt.shape),1),dtype=numpy.float64)
        loadHist[:,0] = lvdt
        # save to attribute
        self.loadHist = loadHist
        return



class RBS(BH):
    """ 
    reduced beam section
    inherit from BH, since same methods will be used.
    """
    #
    # Attributes (Object Initialization)
    #
    def __init__(self, odbPath, material, failureLoad,
                        setName='CenterNode',
                        loadSetName='LVDT'):
        """ override BH init. return instance with desired attributes """
        
        # initialize using superclass but with RBS defaults
        superSpecimen.__init__(self, odbPath=odbPath, instanceName="RBS_", 
                                     setName=setName, material=material,
                                     failureLoad=failureLoad,loadSetName=loadSetName)
        return