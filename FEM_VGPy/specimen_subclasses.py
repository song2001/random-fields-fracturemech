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
    
    Attributes:
        odbPath      = string location of ODB file. Must be of the form:
                       'C:\\temp\\folder\\odbfile.odb'
        material     = string of material type: 'AP50' or 'AP70HP'
        radius       = float of the value for the notch radius
        failureDispl = list or tuple of the failure displacements (test results)
                       for this specimen
        displSetName = string of the name of the set where displacement is applied
        setName      = string of the name of the set of interest
                       (e.g. set to obtain VGI)
    
    Attributes set by self.calcNodalAvgMonoVGI():
        VGI          = numpy array of entire VGI history
                       rows are associated with abaqus frames
                       columns are associated with the nodes in the set
    
    Attributes set by self.calcElemAvgMonoVGI():
        VGI          = same as above, but values are associated with elements
    
    Attributes set by self.determineFailureVGI():
        failureVGI   = numpy array of the VGI's associated with failure displacements
                       row 1 is failureDispl 1, etc.
                       columns are associated with the nodes in the set

    Attributes set by self.fetchVolume():
        elemVol      = numpy array of element volumes in self.setName
        
    Attributes set by self.fetchMesh():
        nodesCoords  = numpy array of the nodal coordinates
        elemConnect  = numpy array of the elemental connectivity
    
    
    """
    
    #
    # Attributes (object initialization)
    #
    def __init__(self, odbPath, material, radius, failureDispl, 
                       setName='CenterNode',
                       displSetName='DisplacementSurface'):
        """ return object with desired attributes """
        
        # initialize the generic superclass variables
        superSpecimen.__init__(self, odbPath=odbPath, instanceName="SNTT_", 
                                     setName=setName, material=material)
        
        # initialize SNTT-specific variables
        self.radius       = radius
        self.failureDispl = failureDispl
        self.displSetName = displSetName.upper() # must be uppercase
        return
    
    #
    # Methods
    #
    def determineFailureVGI(self):
        """ 
        determine which VGI's correspond to failures.
        this is specific to SNTT because failure is specifically 
        related to the observed displacement
        """
        
        # obtain the displacement history of the simulation
        abqDispl = NodalVariable(self.odbPath, 'U', self.displSetName)
        abqDispl.fetchNodalOutput()
        
        # number of frames in the abaqus history
        nframe = abqDispl.resultData.shape[0]
        
        # determine which frames correspond to failure
        failureIndex = []
        for displ in self.failureDispl:
            # for all failure displacements, locate which frame index it corresponds to.
            # this assumes that the first node (column) is equivalent to the
            # other columns, and that the relevant displacement is in the 
            # 2-direction (AKA y-direction).
            
            # we don't want to compare floats directly, so use a percent error
            percent_err = numpy.absolute( (abqDispl.resultData[:,0,1] - displ)/displ )
            
            # lowest percent error => equivalent to failure displ
            frameind = numpy.argmin(percent_err)
            
            # check to see that it is within tolerance
            if percent_err[frameind] < 0.0005:
                # good. save it.
                failureIndex.append(frameind)
            else:
                # provide warning that something went wrong.
                # append nothing. try to continue.
                print "a failure index could not be located!\n"
        
        # therefore, the VGI's at those frames are the failure VGI's... save them!
        self._save_failureVGI(failureIndex)
        return

class CT(superSpecimen):
    """ a compact tension specimen (ASTM E1820)
    
    Attributes:
    
    Attributes set by :
        failureVGI   = array of the VGI's associated with failure displacements
                       row 1 is failureDispl 1, etc.
                       columns are associated with the nodes in the set
        VGI          = array of entire VGI history
                       rows are associated with abaqus frames
                       columns are associated with the nodes in the set

    Methods:
    """
    
    #
    # Attributes (object initialization)
    #
    def __init__(self, odbPath, material, failureJ1, 
                       setName='CrackExtensionPlane', 
                       crackTipSet='CrackTip', crackName='', stepName='Pull'):
        """ return object with desired attributes """
        
        # initialize generic superclass variables
        superSpecimen.__init__(self, odbPath=odbPath, instanceName="CT_", 
                                     setName=setName, material=material)
        
        # initialize CT-specific variables
        self.failureJ1   = failureJ1
        self.stepName    = stepName
        self.crackName   = crackName
        self.crackTipSet = crackName.upper() # must be uppercase
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
        # check if pre-requisits are properly met
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
            percent_err = numpy.absolute( (setCoords[:,0] - target) / target )
            
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
        
    def determineFailureVGI(self):
        """ determine which VGI's are the failure VGI's  """
        
        # obtain the J1 (last contour) history of the simulation
        abqJ1 = CrackVariable(self.odbName, self.stepName, self.crackName)
        abqJ1.getJintegral()
        abqJ1 = abqJ1.resultData[:,-1]
        # number of total frames in the abaqus history
        nframeHist = abqJ1.shape[0]
        
        # find out which frames represent observed failure
        failureIndex = []
        for J1c in self.failureJ1:
            #for each observed J1c, find the nearest abaqus J1 and it's index
            percent_err = numpy.absolute( (abqJ1 - J1c)/J1c )
            
            # the frame we want has the lowest percent error
            index = numpy.argmin(percent_err)
            
            # index is calculated from the J1 (history variable)
            # but VGI is calculated from a field variable.
            # field variables start at frame 0, while history variables
            # start at the first non-zero frame. So, we must add 1 to the index.
            failureIndex.append(int( index + 1 ))
            
        # therefore, the VGI's at those frames are the failure VGI's... save them!
        self._save_failureVGI(failureIndex)
        return

##### None of the below currently works. Just a skeletal framework.
# class BN(CT):
    # """ blunted notch specimen. inherit from CT """


# class BB(SNTT):
    # """ bolt-bearing specimen. inherit from SNTT """
    # def __init__(self, odbPath, material, radius, failureDispl
                        # setName='',
                        # displSetName=('LVDT_top','LVDT_bottom')):
        # """
        # return instance with desired attributes
        # failureDispl should be the ABAQUS displacement... so, account for symmetry
        # """
        
        # # use SNTT init
        # SNTT.__init__(self, odbPath, material, radius, failureDispl,
                       # setName,displSetName)
                       
        # # we dont want SNTT partName
        # self.partName = 'BB_symm'
        # return
        
    # def determineFailureVGI(self):
        # """ 
        # determine which VGI's correspond to failures.
        # this is specific to BB because failure is specifically 
        # related to the observed LVDT (x2) measurements
        # """
        # raise Exception('failure displacement checks not thought through')
        # # obtain the displacement history of the LVDT sets
        # abqDispl = []
        # for i,dset in enumerate(self.displSetName):
            # # for all defined displacement sets
            # abqDispl.append( NodalVariable(self.odbPath, 'U', dset) )
            # abqDispl[i].fetchNodalOutput()
            # # average this history accross nodes, to obtain 1 observation
            # # per frame, per coordinate direction
            # abqDispl[i].avgNodalOutput()
        
        # # obtain difference between LVDT observations
        # abqDispl = numpy.absolute( abqDispl[0].resultData - abqDispl[1].resultData )
        
        # # number of frames in the abaqus history
        # nframe = abqDispl.shape[0]
                
        # # determine which VGI's are the failure VGI's
        # failureVGI = numpy.zeros((len(self.failureDispl),self.VGI.shape[1]),dtype=numpy.float64)
        
        # row = int(0)
        # alreadySaved = []
        # for frame in range(0,nframe):
            # for displ in self.failureDispl:
                # # this assumes that first node (column) is representative of all displacements
                # # and that the displacement is in the y-direction.
                # # the alreadySaved variable is a hack due to abaqus sometimes
                # # saving frame value twice (we dont want the same frame 2x)
                # if (abqDispl[frame,0,1] == displ) and (displ not in alreadySaved):
                    # failureVGI[row,:] = self.VGI[frame,:]
                    # row += 1
                    # alreadySaved.append(displ)
        
        # self.failureVGI = failureVGI
        # return


# class BH(SNTT):
    # """ bolt-hole specimen. inherit from SNTT """
    
    # def __init__(self, odbPath, material, radius, failureDispl
                        # setName=''
                        # displSetName='LVDT'):
        # """
        # return instance with desired attributes
        # failureDispl should be the ABAQUS displacement... so, account for symmetry
        # """
        
        # # use SNTT init
        # SNTT.__init__(self, odbPath, material, radius, failureDispl,
                       # setName,displSetName)
                       
        # # we dont want SNTT partName
        # self.partName = 'BH_symm'
        # return
        
    # def determineFailureVGI(self):
        # """ 
        # determine which VGI's correspond to failures.
        # this is specific to BH because failure is specifically 
        # related to the observed LVDT measurement
        # """
        # raise Exception('failure displacement checks not thought through')
        # # obtain the displacement history of the LVDT set
        # abqDispl = NodalVariable(self.odbPath, 'U', self.displSetName)
        # abqDispl.fetchNodalOutput()
        # # average this history accross nodes, to obtain 1 observation
        # # per frame, per coordinate direction
        # abqDispl.avgNodalOutput()
        
        # # number of frames in the abaqus history
        # nframe = abqDispl.resultData.shape[0]
                
        # # determine which VGI's are the failure VGI's
        # failureVGI = numpy.zeros((len(self.failureDispl),self.VGI.shape[1]),dtype=numpy.float64)
        
        # row = int(0)
        # alreadySaved = []
        # for frame in range(0,nframe):
            # for displ in self.failureDispl:
                # # this assumes that first node (column) is representative of all displacements
                # # and that the displacement is in the y-direction.
                # # the alreadySaved variable is a hack due to abaqus sometimes
                # # saving frame value twice (we dont want the same frame 2x)
                # if (abqDispl.resultData[frame,0,1] == displ) and (displ not in alreadySaved):
                    # failureVGI[row,:] = self.VGI[frame,:]
                    # row += 1
                    # alreadySaved.append(displ)
        
        # self.failureVGI = failureVGI
        # return


# class RBS(BH):
    # """ reduced beam section. inherit from BH """
    # def __init__(self, odbPath, material, radius, failureDispl
                        # setName=''
                        # displSetName='LVDT'):
        # """
        # return instance with desired attributes
        # failureDispl should be the ABAQUS displacement... so, account for symmetry
        # """
        
        # # use SNTT init (like BH)
        # SNTT.__init__(self, odbPath, material, radius, failureDispl,
                       # setName,displSetName)
                       
        # # but we dont want SNTT partName
        # self.partName = 'RBS_'
        # return