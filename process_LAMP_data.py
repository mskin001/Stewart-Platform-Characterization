import os
import csv
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import lampDataFunc

# %% -----------------------------------------------------------------------------------------
file_name = "Prelim Test 001.csv"
folder = "Characterization Data\Results"

sr = 100 # sample rate

plotResponse = True
plotSorted = True
plotDiff = True
plotDirComp = True
plotSpec = True
plotDiffSpec = True
# %% -----------------------------------------------------------------------------------------
dir_PVA_map = np.array([[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23],
               [26, 32, 38, 27, 33, 39, 28, 34, 40, 29, 35, 41, 30, 36, 42, 31, 37, 43]])

full_file = os.path.join(folder, file_name)
with open(full_file, "r") as file:
    csvreader = csv.reader(file, delimiter=',')
    header = next(csvreader)

raw_data = np.loadtxt(full_file, dtype='float', delimiter=',', skiprows=1)
time_col = header.index("H2C time (us)")
raw_time = raw_data[:,time_col]
exp_time = (raw_time - raw_time[0]) / 100

h2c_start = header.index("H2C surge P")
h2c_cols = np.sum(raw_data[:,h2c_start:h2c_start+18], axis=0).nonzero()
h2c_cols = h2c_cols[0] + h2c_start
num_cols = len(h2c_cols)
num_rows = len(raw_data)

c2h_cols = []
for k in range(num_cols):
    map_idx = np.where(dir_PVA_map[0,:] == h2c_cols[k])
    c2h_cols.append(dir_PVA_map[1,map_idx[0][0]])

h2c_data = np.zeros((num_rows,num_cols))
c2h_data = np.zeros((num_rows,num_cols))
for k in range(num_cols):
    h2c_data[:,k] = raw_data[:,h2c_cols[k]]
    c2h_data[:,k] = raw_data[:,c2h_cols[k]]
    if c2h_data[0,k] != 0:
        c2h_data[:,k] = c2h_data[:,k] - c2h_data[0,k]

# %% -----------------------------------------------------------------------------------------
h2c_sort, c2h_sort, diff_sort = lampDataFunc.testDataSort(h2c_data, c2h_data)

h2cPos = h2c_data[:,0::3]
c2hPos = c2h_data[:,0::3]
diff = (h2cPos - c2hPos)
h2c_spec = sp.fft.fftn(h2cPos, axes=0)
c2h_spec = sp.fft.fftn(c2hPos, axes=0)
diff_spec = sp.fft.fftn(diff, axes=0)
N = len(h2c_spec)
n = np.arange(N)
T = N/sr
freq = n/T

dt = exp_time[1] - exp_time[0]
h2cSpec, h2cFreq, c2hSpec, c2hFreq, diffSpec, diffFreq = lampDataFunc.freqDist(h2cPos, c2hPos, dt)
# %% -----------------------------------------------------------------------------------------
units = ["Pos [m]", "Vel [m/s]", "Acc [m/s^2]", 
         "Angle [rad]", "AngVel [rad/s]", "AngAcc [rad/s^2]"]

if plotResponse == True:
    fig, axs = plt.subplots(3)
    for k in range(num_cols):
        axs[k].plot(exp_time,h2c_data[:,k], label="Commanded")
        axs[k].plot(exp_time,c2h_data[:,k], label="Result")
        axs[k].set_ylabel(units[k])
        axs[k].grid(visible=1,which='major',axis='both')
    axs[0].legend()
    plt.xlabel("Time [s]")

if plotSorted == True:
    fig, sortplt = plt.subplots(3)
    for k in range(num_cols):
        sortplt[k].plot(np.arange(0,num_rows,1),h2c_sort[:,k], label="Commanded")
        sortplt[k].plot(np.arange(0,num_rows,1),c2h_sort[:,k], label="Result")
        sortplt[k].set_ylabel(units[k])
        sortplt[k].grid(visible=1,which="major",axis="both")
    sortplt[0].legend()
    plt.xlabel("Index")

if plotDiff == True:
    fig, diffPlt = plt.subplots(3)
    for k in range(num_cols):
        diffPlt[k].plot(h2c_sort[:,k],diff_sort[:,k])
        diffPlt[k].set_xlabel(units[k])
        diffPlt[k].set_ylabel("Residual")
        diffPlt[k].grid(visible=1,which="major",axis="both")
    
if plotDirComp == True:
    plt.figure()
    plt.plot(h2c_sort,c2h_sort)
    plt.legend(["Pos [m]", "Vel [m/s]", "Acc [m/s^2]"])

if plotSpec == True:
    plt.figure()
    plt.stem(freq,np.abs(h2c_spec), "b", markerfmt=" ", basefmt=" ")
    plt.stem(freq,np.abs(c2h_spec), "r", markerfmt=" ", basefmt=" ")
    plt.xlim((0,6))

# if plotDiffSpec == True:
#     plt.figure()
#     plt.plot(h2cFreq,h2cSpec)

print("Program Complete")
plt.show()