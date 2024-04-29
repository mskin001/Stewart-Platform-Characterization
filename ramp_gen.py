import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

testDOF = [0, 0, 1, 0, 0, 0]
testPoints = np.array([[1, 1], [0.25, .75]]) #([amp], [freqs])
phase = [0, 0] # phase angle in degrees
testLength = 20 # [s]
dt = 0.01 # time step
save_test_files = False
save_file_name = "TP016-He-p5-p5-p1-1p2"
#%% -------------------------------------------------------------------------
t = np.arange(0,testLength,dt)
pos = np.zeros((len(t), np.sum(testDOF)))
iter = 0
for k in range(np.sum(testDOF)):
    Amps = np.linspace(testPoints[0,iter], testPoints[0,iter+1], int(testLength/dt))
    freqs = np.linspace(testPoints[1,iter], testPoints[1,iter+1], int(testLength/dt))
    dph = phase[k] * (np.pi / 180)
    pos[:,k] = Amps * np.sin((2*np.pi*freqs*t) + dph)
    iter = iter + 2

vel = np.gradient(pos, axis=0)
acc = np.gradient(vel, axis=0)

freqResp = sp.fft.fftn(pos)
freq = sp.fft.fftfreq(np.size(freqResp), d=dt)
#%% -------------------------------------------------------------------------
#%% Save multisine in .csv files
test_pos = np.zeros((np.size(pos,0),6))
test_vel = np.zeros((np.size(pos,0),6))
test_acc = np.zeros((np.size(pos,0),6))
ind = np.nonzero(testDOF)

for k in range(np.shape(pos)[1]):
    test_pos[:,ind[0][k]] = pos[:,k]
    test_vel[:,ind[0][k]] = vel[:,k]
    test_acc[:,ind[0][k]] = acc[:,k]

if save_test_files:
    pos_file = save_file_name + " pos" + ".csv"
    vel_file = save_file_name + " vel" + ".csv"
    acc_file = save_file_name + " acc" + ".csv"
    time_file = save_file_name + " time" + ".csv"
    np.savetxt(pos_file,test_pos, delimiter=',')
    np.savetxt(vel_file,test_vel, delimiter=',')
    np.savetxt(acc_file,test_acc, delimiter=',')
    np.savetxt(time_file,t, delimiter=',')

#%% -------------------------------------------------------------------------
lUnits = ["Pos", "Vel", "Acc"]
rUnits = ["Angle", "AngVel", "AngAcc"]
DOFs = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

lines = []
if np.sum(testDOF[0:3]) > 0:
    fig1, axl = plt.subplots(3)
    for k in range(np.sum(testDOF[0:3])):
        axl[0].plot(t,pos[:,k])
        axl[0].set_ylabel(lUnits[0])
        axl[0].grid(visible=1,which='major',axis='both')
        axl[1].plot(t,vel[:,k])
        axl[1].set_ylabel(lUnits[1])
        axl[1].grid(visible=1,which='major',axis='both')
        axl[2].plot(t,acc[:,k])
        axl[2].set_ylabel(lUnits[2])
        axl[2].grid(visible=1,which='major',axis='both')
        lines = np.append(lines, DOFs[ind[0][k]])
    axl[0].legend(lines, loc="upper right")
    axl[2].set_xlabel("Time")

if np.sum(testDOF[3:]) > 0:
    fig2, axr = plt.subplots(3)
    lines = []
    for b in range(np.sum(testDOF[3:])):
        axr[0].plot(t,pos[:,b+k+0])
        axr[0].set_ylabel(rUnits[0])
        axr[0].grid(visible=1,which='major',axis='both')
        axr[1].plot(t,vel[:,b+k+0])
        axr[1].set_ylabel(rUnits[1])
        axr[1].grid(visible=1,which='major',axis='both')
        axr[2].plot(t,acc[:,b+k+0])
        axr[2].set_ylabel(rUnits[2])
        axr[2].grid(visible=1,which='major',axis='both')
        lines = np.append(lines, DOFs[ind[0][b+k+0]])
    axr[0].legend(lines, loc="upper right")
    axr[2].set_xlabel("Time")

plt.figure()
plt.stem(freq,np.abs(freqResp), "b", markerfmt=" ", basefmt=" ", linefmt="blue")
plt.xlim((0,2))
plt.xlabel("Frequency [Hz]")
plt.ylabel("(Amplitude)")

plt.show()