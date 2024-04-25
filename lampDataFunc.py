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
    for k in range(np.shape(hostData)[1]):
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
    # dt = int, the time step between data points
    # Outputs:
    # hostSpec = (MxN) array where M is the number of elements in
    #             the spectra and N is the number of columns in 
    #             hostData
    # controlSpec = (MxN) array where M is the number of elements 
    #                in the spectra and N is the number of columns
    #                in controlData
    # diffSpec = (MxN) array where M is the number of elements 
    #                in the spectra and N is the number of columns
    #                in hostData
    # freq = (1xM) array containing the sample frequencies
    diff = np.zeros(hostData.shape)
    for k in range(np.shape(hostData)[1]):
        diff[:,k] = hostData[:,k] - controlData[:,k]
    hostSpec = sp.fft.fftn(hostData, axes=0)
    contSpec = sp.fft.fftn(controlData, axes=0)
    diffSpec = sp.fft.fftn(diff, axes=0)
    freq = sp.fft.fftfreq(np.size(hostSpec[:,0]), d=dt)
        
    return hostSpec, contSpec, diffSpec, freq

def tfestimate(x, y, dirOfInt, dt): #coupled transfer function estimate

    inputArray = np.zeros((len(x[:,1]),6))
    for k in range(6):
        inputArray[:,k] = x[:,dirOfInt]
    fyx, Pyx = sp.signal.csd(y,inputArray, fs=dt, axis=0)
    fxx, Pxx = sp.signal.welch(inputArray, fs=dt, axis=0)
    H = Pyx / Pxx
    ph = np.angle(H, deg=True)

    return H, ph, fyx, fxx