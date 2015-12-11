"""
Vincente Pericoli
UC Davis
10/15/15

"""

#
# imports
#
import numpy
import sys
import os
import myPaths
sys.path.append( myPaths.PyMATLAB() )
import matlab
import matlab.engine

#
# function defs
#
def _convert_dict_dtypes(dictionary):
    """
    takes in dictionary of numpy/python dtypes, and returns 
    dictionary of nearly equivalent MATLAB dtypes.

    type 'list' is currently unsupported if they are mixed dtype.
    (though this could be easily implemented using MATLAB cell,
     it would not work well in the context of FEM_VGPy)
    """
    dict_out = {}
    
    for key in dictionary.keys():
        # walk through all keys, checking the type
        value = dictionary[key]

        # first if-elif ladder (initial checks)
        if value is None:
            #don't add it to the output
            continue
        elif type(value) is str:
            # matlab can automatically convert str to char
            dict_out[key] = value
            continue
        elif type(value) is list:
            # attempt to convert this to a tuple. element
            # dtypes will be taken care of in 2nd ladder
            value = tuple(value)
        
        # second if-elif ladder (convert dtypes)
        if type(value) is tuple:
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

        elif type(value) is dict:
            # use recursion to take care of this
            dict_out[key] = _convert_dict_dtypes(value)
            
        else:
            # undefined. alert user, then save as-is.
            # if there is a problem, matlab engine will
            # throw the proper exceptions
            print "\n!!! undefined type ",
            print type(value),
            print "... saving anyway\n"
            dict_out[key] = value
    
    return dict_out


def VGIPy(inputData, saveKey):
    """ 
    saves a python FEM_VGPy object instance
    (or list/tuple of FEM_VGPy object instances)
    to a MATLAB database
    
    input:
        inputData = single or list/tuple of VGPy objects
        saveKey   = string name to save the MATLAB structure
    """
    # get CWD
    wd = os.getcwd()
    # change to save dir
    os.chdir( myPaths.saveResults() )
    # start matlab engine
    eng = matlab.engine.start_matlab()
    
    if (type(inputData) is list) or (type(inputData) is tuple):
        #if we have a list or tuple
        #then assume the elements are VGIPy classes
        
        #initialize saveData dict which will be exported
        saveData = {}
        
        for instance in inputData:
            #save instance dicts to saveData dict
            #we want the dict to have matlab dtypes
            saveData[instance.name] = _convert_dict_dtypes( instance.__dict__ )
        
        # now, convert dict to struct
        saveStruct = eng.struct(saveData)       
        
        # transport struct to matlab workspace, and save
        eng.workspace[saveKey] = saveStruct
        eng.save(saveKey + '.mat',saveKey,nargout=0)
        
    else:
        #assume input is a VGIPy class or otherwise has a name attribute
        saveData = {}
        
        try:
            #this works if inputData is an instance
            saveData[inputData.name] = _convert_dict_dtypes( inputData.__dict__ )
        except:
            #this works if inputData is an instance.__dict__
            saveData[inputData["name"]] = _convert_dict_dtypes( inputData )

        # now, convert dict to struct
        saveStruct = eng.struct(saveData)       
        
        # transport struct to matlab workspace, and save
        eng.workspace[saveKey] = saveStruct
        eng.save(saveKey + '.mat',saveKey,nargout=0)


    # exit matlab engine
    eng.quit()    
    # return to previous dir
    os.chdir( wd )
    # alert user
    print "MATLAB Binary Database saved to: " + myPaths.saveResults()
    return

