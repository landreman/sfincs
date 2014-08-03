def readHDF5(variableName):
    # First, check whether sfincsOutput.h5 exists:
    
    import os.path
    
    if not os.path.isfile("sfincsOutput.h5"):
        print "Error! The file sfincsOutput.h5 has not been created."
        exit(1)
        
    # Try executing h5dump with no arguments, just to see if h5dump is installed:
        
    import subprocess
        
    try:
        p = subprocess.Popen("h5dump", stdout=subprocess.PIPE)
    except:
        print "Error! Unable to execute h5dump."
        exit(1)

    # If we made it this far, then the h5dump command is indeed accessible to the shell. Next call h5dump with arguments:

    try:
        p = subprocess.Popen(["h5dump", "-y", "-d", "/run  1/"+variableName, "sfincsOutput.h5"], stdout=subprocess.PIPE)
    except:
        print "Error! Unable to read sfincsOutput.h5 using h5dump. It is likely that sfincsOutput.h5 is corrupted."
        exit(1)

    (output, err) = p.communicate()

    # Pick out the line with the variable contents:
    resultsString = output.splitlines()[5]
    # Convert to an array of strings, then to an array of floats:
    return map(float,resultsString.split(","))

def shouldBe(variableName, index, trueValue, relativeTolerance):
    latestValue = readHDF5(variableName)[index]
    relativeDifference = abs((latestValue - trueValue) / trueValue)
    if relativeDifference > relativeTolerance:
        print "*** TEST FAILED!!  Variable "+variableName+"["+str(index)+"] should be close to "+str(trueValue)+", but it is instead "+str(latestValue)
        return 1
    else:
        print "    Test passed:   Variable "+variableName+"["+str(index)+"] should be close to "+str(trueValue)+", and it came out to be "+str(latestValue)+", which is within tolerance."
        return 0
  
