import numpy as np
from matplotlib import pyplot as plt
import scipy as sp

#Amp = [1] # maximum amplitude in [m]
Amp = [0.1, 1.5]
#freq = [0.1, 2] # test freqency range [minFreq, maxFeq]
freq = [0.1]
phase = 0 # phase angle in degrees
testLength = 10 # [s]
dt = 0.01 # time step

#%% -------------------------------------------------------------------------
t = np.arange(0,testLength,dt)

if len(Amp) > 1 and len(freq) > 1:
    print("Max amplutde and frequency cannot both vary")
    exit()
elif len(Amp) == 1 and len(freq) > 1:
    freqs = np.linspace(freq[0], freq[1], int(testLength/dt)) # construct frequency array
    dph = phase * (np.pi / 180) # convert to rad
    pos = Amp * np.sin((2*np.pi*freqs*t) + dph)
elif len(Amp) > 1 and len(freq) == 1:
    Amps = np.linspace(Amp[0], Amp[1], int(testLength/dt)) # construct frequency array
    dph = phase * (np.pi / 180) # convert to rad
    pos = np.linspace(0, 1, len(t)) * np.sin((2*np.pi*freq[0]*t) + dph)

vel = np.gradient(pos)
acc = np.gradient(vel)

posFreq = sp.fft.fft(pos)
freqBins = sp.fft.fftfreq(np.size(posFreq), d=dt)

#%% -------------------------------------------------------------------------
fig1, ax1 = plt.subplots(3)
ax1[0].plot(t,pos)
ax1[0].set_ylabel("Pos")
ax1[0].grid(visible=1,which='major',axis='both')
ax1[1].plot(t,vel)
ax1[1].set_ylabel("Vel")
ax1[1].grid(visible=1,which='major',axis='both')
ax1[2].plot(t,acc)
ax1[2].set_ylabel("Acc")
ax1[2].set_xlabel("Time [s]")
ax1[2].grid(visible=1,which='major',axis='both')

# plt.figure()
# plt.plot(t,freqs)

plt.figure()
plt.stem(freqBins,np.abs(posFreq), basefmt=" ", markerfmt=" ")
plt.xlim(0,5)

plt.show()