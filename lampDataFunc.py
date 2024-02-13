# LAMP Data processing functions. This file was created to house any
# custom functions created for the LAMP characterization data
# processing. For notation, host data is what we send to the
# platform. This is the demanded motion, i.e. the signal produced by
# signal_gen.py. The control data is the responce from the platform.
# This is what we get back as actual motion.

import numpy as np
import scipy as sp

def testDataSort(hostData, controlData):
    # Sorts the host and controller data from smallest to largest. 
    # This function also finds the difference (i.e. residual)
    # between the host and controller data. It returns the sorted
    # data and the resitual for each direction. 
    # Inputs:
    # hostData = (MxN) array where M is the length of the test
    #             and N is the 3x the number of DOFs in the test
    # controlData = (MxN) array where M is the length of the test
    #                 and N is the 3x the number of DOFs in the 
    #                 test
    #
    # Outputs:
    # hostSort = (MxN) sortedarray where M is the length of the 
    #             test and N is the 3x the number of DOFs in the
    #             test.
    hostSort = np.zeros(hostData.shape)
    contSort = np.zeros(controlData.shape)
    diff = np.zeros(hostData.shape)
    numCols = np.shape(hostData)[1]
    for k in range(numCols):
        hostSort[:,k] = np.sort(hostData[:,k])
        contSort[:,k] = np.sort(controlData[:,k])
        diff[:,k] = contSort[:,k] - hostSort[:,k]
    return hostSort, contSort, diff

def freqDist(hostData, controlData):
    # Computes the fft column wise on the host and controller data
    # and returns the frequency spectra.
    # Inputs:
    # hostData = (MxN) array where M is the length of the test
    #             and N is the 3x the number of DOFs in the test
    # controlData = (MxN) array where M is the length of the test
    #                 and N is the 3x the number of DOFs in the 
    #                 test
    # Outputs:
    # hostSpec = (MxN) array where M is the number of elements in
    #             the spectra and N is the number of columns in 
    #             hostData
    # controlSpec = (MxN) array where M is the number of elements 
    #                in the spectra and N is the number of columns
    #                in hostData
    hostSpec = []
    contSpec = []
    for k in range(np.shape(hostData)[1]):
        hostSpec = sp.fft.fft(hostData[:,k])
        contSpec = sp.fft.fft(controlData[:,k])
    return hostSpec, contSpec


