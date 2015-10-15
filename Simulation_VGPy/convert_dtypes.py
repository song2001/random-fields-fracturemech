"""
Vincente Pericoli
UC Davis
10/12/15

Using MATLAB PyEngine, convert dtypes into MATLAB equivalents
"""
import numpy
import sys
import myPaths
sys.path.append( myPaths.PyMATLAB() )
import matlab

def convert_dict_dtypes(dictionary):
    """
    takes in dictionary of numpy types, and returns dictionary
    of nearly equivalent MATLAB types.

    type 'list' is currently unsupported
    (list items can be of arbitrary type... need a matlab cell or similar)
    """
    dict_out = {}
    
    for key in dictionary.keys():
        #walk through all keys, checking the type
        value = dictionary[key]

        if value is None:
            #don't add it to the output
            pass

        elif type(value) is str:
            # matlab can automatically convert str to char
            pass
        
        elif type(value) is tuple:
            #iterable type; convert to array
            if type(value[0]) is int:
                dict_out[key] = matlab.int32(value)
            elif type(value[0]) is float:
                dict_out[key] = matlab.double(value)

        elif type(value) is numpy.ndarray:
            # we want to convert to an equivalent matlab matrix...

            #convert data to equivalent list
            data_list = value.tolist()
            
            #change dtype to matlab double and save
            dict_out[key] = matlab.double(data_list)
            
        elif type(value) is float:
            # matlab considers floats to be matrices as well.
            dict_out[key] = matlab.double([value])

        elif type(value) is int:
            # matlab considers ints to be matrices as well.
            dict_out[key] = matlab.int32([value])

        else:
            # undefined. alert user, then save as-is.
            # if there is a problem, matlab engine will
            # throw the proper exceptions
            print "\n!!! undefined type ",
            print type(value),
            print "... saving anyway\n"
            dict_out[key] = value
    
    return dict_out
