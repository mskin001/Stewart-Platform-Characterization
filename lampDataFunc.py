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

def freqDist(hostData, controlData, dt):
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
    hostSpec = sp.fft.fftn(hostData, axes=0)
    contSpec = sp.fft.fftn(controlData, axes=0)
    diff = hostData - controlData
    diffSpec = sp.fft.fftn(diff, axis=0)

    hostFreq = np.zeros(np.shape(hostSpec))
    contFreq = np.zeros(np.shape(contSpec))
    diffFreq = np.zeros(np.shape(diffSpec))
    for k in range(np.shape(hostSpec)[1]):
        hostFreq[:,k] = sp.fft.fftfreq(hostSpec[:,k], d=dt)
        contFreq[:,k] = sp.fft.fftfreq(contSpec[:,k], d=dt)
        diffFreq[:,k] = sp.fft.fftfreq(diffSpec[:,k], d=dt)
        
    return hostSpec, hostFreq, contSpec, contFreq, diffSpec, diffFreq


