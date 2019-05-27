#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main regcoil directory.

execfile('../testsCommon.py')

numFailures = 0

outputFile = readOutputFile()
referenceFile = readReferenceFile()

absoluteTolerance = 1e-13
relativeTolerance = 1e-100
numFailures += compareToReference(referenceFile,outputFile,'area_plasma',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'area_coil',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'volume_plasma',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'volume_coil',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'ntheta_plasma',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'nzeta_coil',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'ntheta_plasma',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'nzeta_coil',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'mpol_magnetization',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'ntor_magnetization',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'num_basis_functions',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'system_size',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'lambda',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'s_magnetization',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'s_integration',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'s_weights',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'RHS_B',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'mean_curvature_coil',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'Jacobian_coil',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance)
numFailures += compareToReference(referenceFile,outputFile,'chi2_M',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance,skipFirstElement=True)
numFailures += compareToReference(referenceFile,outputFile,'chi2_B',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance,skipFirstElement=True)
numFailures += compareToReference(referenceFile,outputFile,'max_M',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance,skipFirstElement=True)
numFailures += compareToReference(referenceFile,outputFile,'min_M',absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance,skipFirstElement=True)

outputFile.close()
referenceFile.close()

print "numFailures:",numFailures
exit(numFailures > 0)
