import os
import csv
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import lampDataFunc

# %% -----------------------------------------------------------------------------------------
#file_name = "RW_Second_Test_UDP.csv"
file_name = "EM001_010-6_2.csv"
emfolder = "Characterization Data\Emulator Results"
tpfolder = "Characterization Data\Test Profiles"
rwfolder = "Characterization Data\Real World Results"

sr = 100 # sample rate
dt = 0.01
DOF = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

plotResponse = True
plotDiff = True
plotSorted = True
plotSortedDiff = True
plotDirComp = False
plotSpec = True
plotDiffSpec = False
# %% -----------------------------------------------------------------------------------------
dir_PVA_map = np.array([[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
               [26, 32, 38, 27, 33, 39, 28, 34, 40, 29, 35, 41, 30, 36, 42, 31, 37, 43]])

if file_name[0:2] == "EM":
    folder = emfolder
elif file_name[0:2] == "RW":
    folder = rwfolder
else:
    folder = tpfolder

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
diff = np.zeros(h2c_data.shape)
for k in range(num_cols):
    diff[:,k] = (h2c_data[:,k] - c2h_data[:,k])

h2c_sort, c2h_sort, diff_sort = lampDataFunc.testDataSort(h2c_data, c2h_data)

h2cPos = h2c_data[:,0::3]
c2hPos = c2h_data[:,0::3]
h2cSpec, c2hSpec, diffSpec, freq = lampDataFunc.freqDist(h2cPos, c2hPos, dt)
# %% -----------------------------------------------------------------------------------------
units = ["Pos [m]", "Vel [m/s]", "Acc [m/s^2]", 
         "Angle [rad]", "AngVel [rad/s]", "AngAcc [rad/s^2]"]

if plotResponse == True:
    iter = 0
    for b in range(int(num_cols/3)):
        fig, axs = plt.subplots(3)
        for k in range(3):
            col = k + iter
            axs[k].plot(exp_time,h2c_data[:,col], label="Commanded")
            axs[k].plot(exp_time,c2h_data[:,col], label="Result")
            axs[k].set_ylabel(units[k])
            axs[k].grid(visible=1,which='major',axis='both')
        iter = iter + 3
        axs[0].set_title(DOF[b])
        axs[0].legend()
    plt.xlabel("Time [s]")

if plotDiff == True:
    iter = 0
    for b in range(int(num_cols/3)):
        fig, diffPlt = plt.subplots(3)
        for k in range(3):
            col = k + iter
            diffPlt[k].plot(exp_time,diff[:,col], label="Commanded")
            diffPlt[k].set_ylabel(units[k])
            diffPlt[k].grid(visible=1,which='major',axis='both')
        iter = iter + 3
        axs[0].set_title(DOF[b])
    plt.xlabel("Time [s]")

if plotSorted == True:
    iter = 0
    for b in range(int(num_cols/3)):
        fig, sortplt = plt.subplots(3)
        for k in range(3):
            col = k + iter
            sortplt[k].plot(np.arange(0,num_rows,1),h2c_sort[:,col], label="Commanded")
            sortplt[k].plot(np.arange(0,num_rows,1),c2h_sort[:,col], label="Result")
            sortplt[k].set_ylabel(units[k])
            sortplt[k].grid(visible=1,which="major",axis="both")
        iter = iter + 3
        sortplt[0].set_title(DOF[b])
        sortplt[0].legend()
    plt.xlabel("Index")

if plotSortedDiff == True:
    fig, sortDiff = plt.subplots(3)
    for k in range(num_cols):
        sortDiff[k].plot(h2c_sort[:,k],diff_sort[:,k])
        sortDiff[k].set_xlabel(units[k])
        sortDiff[k].set_ylabel("Residual")
        sortDiff[k].grid(visible=1,which="major",axis="both")
    
if plotDirComp == True:
    plt.figure()
    plt.plot(h2c_sort,c2h_sort)
    plt.legend(["Pos [m]", "Vel [m/s]", "Acc [m/s^2]"])

if plotSpec == True:
    for k in range(np.shape(h2cSpec)[1]):
        plt.figure()
        plt.stem(freq,np.abs(h2cSpec[:,k]), "b", markerfmt=" ", basefmt=" ")
        plt.stem(freq,np.abs(c2hSpec[:,k]), "r", markerfmt=" ", basefmt=" ")
        plt.title(DOF[k])
        plt.xlim((0,6))
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("(Amplitude)")

if plotDiffSpec == True:
    plt.figure()
    plt.stem(freq,np.abs(diffSpec), markerfmt=" ", basefmt=" ")
    plt.xlim(0,6)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("(Amplitude)")

print("Program Complete")
plt.show()