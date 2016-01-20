"""
Vincente Pericoli
UC Davis

Calculates the Void Growth Index for micromechanical fracture
"""

# imports
import numpy

# function definitions
def __check_input_args(mises, pressure, PEEQ):
    """ performs a rudimentary check on inputs """
    
    # rudimentary check on inputs
    if not ((mises.shape == pressure.shape) and (pressure.shape == PEEQ.shape)):
        raise Exception('calcVGI: Matrices are not the same shape!')
    elif not ((type(mises) == type(pressure)) and (type(pressure) == type(PEEQ))):
        raise Exception('calcVGI: inputs are not the same type!')
    elif not isinstance(mises, numpy.ndarray):
        raise Exception('calcVGI: inputs are not of type numpy.ndarray!')      
    return

def calcMonotonicVGI(mises, pressure, PEEQ):
    """
    Takes rank-2 matrices of mises, pressure, PEEQ
    returns a rank-2 matrix of monotonic VGI
    
    matrix must be ordered such that the rows are different
    "history" values, and the columns are different "nodes"
    or other such distinctly different objects
    
    Verified to produce accurate results: 09/21/2015
    """
    
    # check input args
    __check_input_args(mises, pressure, PEEQ)
    
    # preallocate triax
    triax = numpy.zeros(mises.shape, dtype=numpy.float64)
    
    # calculate the stress triaxiality
    # stress triaxiality is an element-wise divide of -pressure/mises
    # skip first element (since division by zero), this will remain 0
    triax[1:,:] = -pressure[1:,:]/mises[1:,:]
    
    # calculate the integrand
    integrand = numpy.exp(1.5*triax)
    
    # preallocate
    VGI  = numpy.zeros(mises.shape, dtype=numpy.float64)
    nrow = mises.shape[0]
    ncol = mises.shape[1]
    
    # trap rule numerical integration:
    for row in range(1,nrow):
        # for all rows ("history" values), except first row
        for col in range(0,ncol):
            # for all columns (or "objects" with distinctly different VGIs)
            
            # incremental VGI
            dVGI = 0.5*(PEEQ[row,col] - PEEQ[row-1,col])*(integrand[row,col] + integrand[row-1,col])
            
            # sum into VGI
            VGI[row,col] = VGI[row-1,col] + dVGI

    return VGI

    
def calcCyclicVGI(mises, pressure, PEEQ):
    """
    Input: rank-2 matrices of mises, pressure, PEEQ.
    
    Output: a tuple of
    (rank-2 matrix of VGI, rank-2 matrix of cumulative damage PEEQ)
    
    input matrices must be ordered such that the rows are different
    "history" values, and the columns are different "nodes"
    or other such distinctly different objects.
    
    This also works fine for monotonic loading, though it would
    perform unnecessary calcs and the cumulative damage PEEQ 
    output is meaningless in that context.
    """
    
    # check input args
    __check_input_args(mises, pressure, PEEQ)
    
    # determine problem size
    nrow = mises.shape[0]
    ncol = mises.shape[1]
    
    # preallocate arrays
    triax    = numpy.zeros(mises.shape, dtype=numpy.float64)
    VGI      = numpy.zeros(mises.shape, dtype=numpy.float64)
    cumePEEQ = numpy.zeros(mises.shape, dtype=numpy.float64)
    
    # calculate the stress triaxiality
    # stress triaxiality is an element-wise divide of -pressure/mises
    # skip first element (since division by zero), this will remain 0
    triax[1:,:] = -pressure[1:,:]/mises[1:,:]
    
    # calculate the integrand (with absolute value of triax)
    integrand = numpy.exp(1.5*numpy.absolute(triax))
    
    # calculate VGI and damage
    for row in range(1,nrow):
        # for all rows ("history" values), except first row (initial zero-frame)
        for col in range(0,ncol):
            # for all columns (or "objects" with distinctly different VGIs)
            
            # calculate incremental VGI (trap rule numerical integration)
            dPEEQ = PEEQ[row,col] - PEEQ[row-1,col]
            dVGI = 0.5 * dPEEQ * (integrand[row,col] + integrand[row-1,col])
            dVGI = dVGI * numpy.sign( triax[row,col] )
            
            # sum into VGI
            VGI[row,col] = VGI[row-1,col] + dVGI
            if VGI[row,col] <= 0.0:
                # VGI can't be <= 0
                VGI[row,col] = 0.0
            
            # calculate the corresponding damage
            if triax[row,col] < 0.0:
                # compressive excursion, damage occurs
                cumePEEQ[row,col] = cumePEEQ[row-1,col] + dPEEQ
            else:
                # tensile excursion, normal VGI behavior
                cumePEEQ[row,col] = cumePEEQ[row-1,col]

    raise Exception('Damage cycle-counting has not been checked with Myers et al 2009')
    return (VGI, cumePEEQ)

