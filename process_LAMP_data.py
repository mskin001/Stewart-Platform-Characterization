import os
import csv
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
import lampDataFunc

# %% -----------------------------------------------------------------------------------------
#file_name = "Pioneer_2024_02_13_20_01_38.csv"
file_name = "EM001_010-6_2.csv"
emfolder = "Characterization Data\Emulator Results"
tpfolder = "Characterization Data\Test Profiles"
rwfolder = "Characterization Data\Real World Results"

sf = 25 # sample rate
dt = 0.01
DOF = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

plotResponse = True
plotDiff = True
plotSorted = False
plotSortedDiff = False
plotDirComp = False
plotSpec = False
plotDiffSpec = False
plotBode = True
# %% -----------------------------------------------------------------------------------------
dir_PVA_map = np.array([[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
               [26, 32, 38, 27, 33, 39, 28, 34, 40, 29, 35, 41, 30, 36, 42, 31, 37, 43]])

if file_name[0:2] == "EM":
    folder = emfolder
elif file_name[0:2] == "RW":
    folder = rwfolder
else:
    #folder = tpfolder
    folder = "C:\\Users\\mskinner\\NREL\\Water Power Equipment - Motion Platform (LAMP)\\LAMP Characterization"

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

wind = ("boxcar", 2)
tf, ph, fyx, fxx = lampDataFunc.tfestimate(h2cPos, c2hPos, len(h2cPos)/sf)

# %% -----------------------------------------------------------------------------------------
units = ["Pos [m]", "Vel [m/s]", "Acc [m/s^2]", 
         "Angle [rad]", "AngVel [rad/s]", "AngAcc [rad/s^2]"]

if plotResponse == True:
    iter = 0
    for b in range(int(num_cols/3)):
        fig, axs = plt.subplots(3, sharex=True)
        for k in range(3):
            col = k + iter
            axs[k].plot(exp_time,h2c_data[:,col], label="Input")
            axs[k].plot(exp_time,c2h_data[:,col], label="Output")
            axs[k].set_ylabel(units[k])
            axs[k].grid(visible=1,which='major',axis='both')
        iter = iter + 3
        axs[0].set_title(DOF[b])
        axs[0].legend()
    plt.xlabel("Time [s]")

if plotDiff == True:
    iter = 0
    for b in range(int(num_cols/3)):
        fig, diffPlt = plt.subplots(3, sharex=True)
        for k in range(3):
            col = k + iter
            diffPlt[k].plot(exp_time,diff[:,col], label="Residual")
            diffPlt[k].set_ylabel(units[k])
            diffPlt[k].grid(visible=1,which='major',axis='both')
        iter = iter + 3
        axs[0].set_title(DOF[b])
    plt.xlabel("Time [s]")

if plotSorted == True:
    iter = 0
    for b in range(int(num_cols/3)):
        fig, sortplt = plt.subplots(3, sharex=True)
        for k in range(3):
            col = k + iter
            sortplt[k].plot(np.arange(0,num_rows,1),h2c_sort[:,col], label="Input")
            sortplt[k].plot(np.arange(0,num_rows,1),c2h_sort[:,col], label="Output")
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
        plt.stem(freq*10,np.abs(h2cSpec[:,k]), "b", markerfmt=" ", basefmt=" ", linefmt="blue")
        plt.stem(freq*10,np.abs(c2hSpec[:,k]), "r", markerfmt=" ", basefmt=" ", linefmt="orange")
        plt.title(DOF[k])
        plt.xlim((0,0.5))
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("(Amplitude)")

if plotDiffSpec == True:
    plt.figure()
    plt.stem(freq*10,np.abs(diffSpec), markerfmt=" ", basefmt=" ")
    plt.xlim(0,0.5)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("(Amplitude)")

if plotBode == True:
    fig, bod = plt.subplots(2, sharex=True)
    bod[0].semilogx(fyx*10, (np.abs(tf.T)))
    bod[1].semilogx(fyx*10, ph.T)
    bod[0].set_ylabel("Gain")
    bod[0].grid(visible=1,which="major",axis="both")
    bod[1].grid(visible=1,which="major",axis="both")
    bod[1].set_ylabel("Phase [deg]")
    plt.xlabel("Frequency [Hz]")

print("Program Complete")
plt.show()