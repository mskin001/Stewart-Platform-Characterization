import numpy as np
import scipy as sp

def testDataSort(hostData, controlData):
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
    hostSpec = []
    contSpec = []
    for k in range(np.shape(hostData)[1]):
        hostSpec = sp.fft.fft(hostData[:,k])
        contSpec = sp.fft.fft(controlData[:,k])
    return hostSpec, contSpec


