import os
import csv
from matplotlib import pyplot as plt
import numpy as np
import lampDataFunc

# %% -----------------------------------------------------------------------------------------
file_name = "Prelipipm Test 001.csv"
folder = "C:\\Users\skinn\Documents\GitHub\Stewart-Platform-Characterization\Characterization Data\Results"
#folder = "C:\\Users\\mskinner\\NREL\Water Power Equipment - Motion Platform (LAMP)\LAMP Characterization\Characterization Data\Results"

plotResponse = True
plotSorted = True
plotDiff = True
plotDirComp = True
# %% -----------------------------------------------------------------------------------------
dir_PVA_map = np.array([[4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 18, 19, 20, 21, 22, 23],
               [26, 32, 38, 27, 33, 39, 28, 34, 40, 29, 35, 41, 30, 36, 42, 31, 37, 43]])

full_file = os.path.join(folder, file_name)
print(" ")
print(full_file)
print(" ")
with open(full_file, "r") as file:
    csvreader = csv.reader(file, delimiter=',')
    header = next(csvreader)

raw_data = np.loadtxt(full_file, dtype='float', delimiter=',', skiprows=1)
time_col = header.index("H2C time (us)")
exp_time = raw_data[:,time_col]

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

# %% -----------------------------------------------------------------------------------------
h2c_sort, c2h_sort, diff = lampDataFunc.testDataSort(h2c_data, c2h_data)

# %% -----------------------------------------------------------------------------------------
if plotResponse == True:
    fig, axs = plt.subplots(3)
    for k in range(num_cols):
        axs[k].plot(exp_time,h2c_data[:,k], label="Commanded")
        axs[k].plot(exp_time,c2h_data[:,k], label="Result")
    axs[0].set_ylabel("Amplitude [m]")
    axs[0].grid(visible=1,which='major',axis='both')
    axs[1].set_ylabel("Vel [m/s]")
    axs[1].grid(visible=1,which='major',axis='both')
    axs[2].set_ylabel("Acc [m/s^2]")
    axs[2].grid(visible=1,which='major',axis='both')
    plt.xlabel("Time [s]")
    axs[0].legend()

if plotSorted == True:
    fig, sortplt = plt.subplots(3)
    for k in range(num_cols):
        sortplt[k].plot(np.arange(0,num_rows,1),h2c_sort[:,k], label="Commanded")
        sortplt[k].plot(np.arange(0,num_rows,1),c2h_sort[:,k], label="Result")
    sortplt[0].set_ylabel("Position [m]")
    sortplt[0].grid(visible=1,which="major",axis="both")
    sortplt[0].legend()
    sortplt[1].set_ylabel("Velocity [m/s]")
    sortplt[1].grid(visible=1,which="major",axis="both")
    sortplt[2].set_ylabel("Acceleration [m/s^2]")
    sortplt[2].grid(visible=1,which="major",axis="both")
    plt.xlabel("Index")

if plotDiff == True:
    fig, diffPlt = plt.subplots(3)
    for k in range(num_cols):
        diffPlt[k].plot(h2c_sort[:,k],diff[:,k])
    diffPlt[0].set_xlabel("Pos [m]")
    diffPlt[0].set_ylabel("Residual")
    diffPlt[0].grid(visible=1,which="major",axis="both")
    diffPlt[1].set_xlabel("Vel [m/s]")
    diffPlt[1].set_ylabel("Residual")
    diffPlt[1].grid(visible=1,which="major",axis="both")
    diffPlt[2].set_xlabel("Acc [m/s^2]")
    diffPlt[2].set_ylabel("Residual")
    diffPlt[2].grid(visible=1,which="major",axis="both")
    
if plotDirComp == True:
    plt.figure()
    plt.plot(h2c_sort,c2h_sort)
    plt.legend(["Pos [m]", "Vel [m/s]", "Acc [m/s^2]"])

print("Program Complete")
plt.show()