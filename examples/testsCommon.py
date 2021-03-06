import os
import numpy as np

def readOutputFile():

    head, dirname = os.path.split(os.getcwd())
    outputFilename = "regcoil_out."+dirname+".nc"
    
    if not os.path.isfile(outputFilename):
        print "Error! The output file "+outputFilename+" has not been created."
        exit(1)
        
    from scipy.io import netcdf
    try:
        f = netcdf.netcdf_file(outputFilename,'r')
    except:
        print "ERROR! Unable to read netCDF output file "+outputFilename
        raise

    print "Reading test output file "+outputFilename
    return f

def readReferenceFile():

    head, dirname = os.path.split(os.getcwd())
    outputFilename = "regcoil_out."+dirname+".reference.nc"
    
    if not os.path.isfile(outputFilename):
        print "Error! The reference output file "+outputFilename+" cannot be found."
        exit(1)
        
    from scipy.io import netcdf
    try:
        f = netcdf.netcdf_file(outputFilename,'r')
    except:
        print "ERROR! Unable to read netCDF reference output file "+outputFilename
        raise

    print "Reading reference output file "+outputFilename
    return f

def shouldBe(latestValue, trueValue, relativeTolerance, absoluteTolerance):
    difference = abs(latestValue-trueValue)
    relativeDifference = 0
    if abs(trueValue) > 0:
        relativeDifference = abs(difference / trueValue)
        relativeTest = (relativeDifference <= relativeTolerance)
    else:
        relativeTest = False
    absoluteTest = (difference <= absoluteTolerance)
    #string = "Expected a value close to "+str(trueValue)+", and it was "+str(latestValue)+". Abs diff="+str(difference)+", rel diff="+str(relativeDifference)
    string = "Expected a value close to {:22.15e}, and it was {:22.15e}. Abs diff={:22.15e}, rel diff={:22.15e}".format(trueValue,latestValue,difference,relativeDifference)
    if relativeTest:
        if absoluteTest:
            print "    Test passed. "+string+". Both abs and rel tol met."
            return 0
        else:
            print "    Test passed. "+string+". Rel tol met. Abs tol not met."
            return 0
    else:
        if absoluteTest:
            print "    Test passed. "+string+". Abs tol met. Rel tol not met."
            return 0
        else:
            print "*** TEST FAILED! "+string+". Neither rel nor abs tol met."
            return 1


def arrayShouldBe(variableName,latestValues, trueValues, relativeTolerance, absoluteTolerance, requireSameLength = True, skipFirstElement = False):
    print "  Comparing "+variableName
    # These next few lines are a hack so this function can be called on scalars without an exception
    try:
        temp = len(latestValues)
    except:
        latestValues = np.array([latestValues])

    try:
        temp = len(trueValues)
    except:
        trueValues = np.array([trueValues])

    if skipFirstElement:
        latestValues = latestValues[1:]
        trueValues = trueValues[1:]

    if requireSameLength and (len(latestValues) != len(trueValues)):
        print "*** TEST FAILED!! Variable "+variableName+" should have length "+str(len(trueValues))+" but it instead has length "+str(len(latestValues))
        return 1

    if len(latestValues) < len(trueValues):
        print "*** TEST FAILED!! Variable "+variableName+" should have length at least "+str(len(trueValues))+" but it instead has length "+str(len(latestValues))
        return 1

    #latestValuesFlat = np.flatten(latestValues)
    #trueValuesFlat = np.flatten(trueValues)
    latestValuesFlat = latestValues.flatten()
    trueValuesFlat = trueValues.flatten()

    numArrayErrors = 0
    for i in range(len(trueValuesFlat)):
        numArrayErrors += shouldBe(latestValuesFlat[i],trueValuesFlat[i],relativeTolerance, absoluteTolerance)

    return numArrayErrors

def compareToReference(referenceFile,testFile,varName,relativeTolerance=1.0e-100,absoluteTolerance=1.0e-13,skipFirstElement=False):
    return arrayShouldBe(varName,testFile.variables[varName][()], referenceFile.variables[varName][()], relativeTolerance, absoluteTolerance, skipFirstElement=skipFirstElement)
