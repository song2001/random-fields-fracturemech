"""
Vincente Pericoli
UC Davis

Calculates the Void Growth Index for micromechanical fracture
Verified to produce accurate results: 09/21/2015
"""

import numpy

def calcMonotonicVGI(mises, pressure, PEEQ):
    """
    Takes rank-2 matrices of mises, pressure, PEEQ
    returns a rank-2 matrix of monotonic VGI
    
    matrix must be ordered such that the rows are different
    "history" values, and the columns are different "nodes"
    or other such distinctly different objects
    """

    #rudimentary check on inputs
    if not ((mises.shape == pressure.shape) and (pressure.shape == PEEQ.shape)):
        raise Exception('calcMonotonicVGI: Matrices are not the same shape!')
    elif not ((type(mises) == type(pressure)) and (type(pressure) == type(PEEQ))):
        raise Exception('calcMonotonicVGI: inputs are not the same type!')
    elif not isinstance(mises, numpy.ndarray):
        raise Exception('calcMonotonicVGI: inputs are not of type numpy.ndarray!')

    # preallocate triax
    triax = numpy.zeros(mises.shape, dtype=numpy.float64)
    
    # calculate the stress triaxiality
    # stress triaxiality is an element-wise divide of -pressure/mises
    triax[1:-1,:] = -pressure[1:-1,:]/mises[1:-1,:]
    triax[0,:] = 0.0 #set to avoid NaN (since division by zero)
    
    # calculate the integrand
    integrand = numpy.exp(1.5*triax)
    
    # preallocate
    VGI  = numpy.zeros(mises.shape, dtype=numpy.float64)
    nrow = mises.shape[0]
    ncol = mises.shape[1]
    
    #trap rule numerical integration:
    for row in range(1,nrow):
        #for all rows, except first row
        for col in range(0,ncol):
            #for all columns
            
            #incremental VGI
            dVGI = 0.5*(PEEQ[row,col] - PEEQ[row-1,col])*(integrand[row,col] + integrand[row-1,col])
            
            #sum into VGI
            VGI[row,col] = VGI[row-1,col] + dVGI

    return VGI
