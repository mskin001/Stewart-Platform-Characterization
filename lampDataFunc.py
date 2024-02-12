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
        hostTemp = sp.fft.fft(hostData[:,0])
        np.append(hostSpec,hostTemp,axis=0)
        contTemp = sp.fft.fft(controlData[:,0])
        np.append(contSpec,contTemp,axis=0)
    print(np.shape(hostTemp))
    print(np.shape(contTemp))
    return hostSpec, contSpec


