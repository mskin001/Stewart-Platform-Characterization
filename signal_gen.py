import numpy as np
import itertools as itt
from matplotlib import pyplot as plt
import scipy as sp

#%% Initialize parameters
test_DOF = np.array([1, 0, 1, 0, 0, 0]) #[surge, sway, heave, roll, pitch, yaw]
save_file_name = "Prelim 007"
save_test_files = False

T = 180
df = 1/T
dt = 0.01
f_range = [0.1, 6]
c = 2 # exponential factor controlling random noise decay, 1 = pink noise
peak_ampl = 2 #approximate peak amplitude, subject to change based on randomness
atten = 0.75
numPhases = 100
reps = 1

#%% Define frequency, omega, and time vectors
numDOF = sum(test_DOF)

f_min = np.round(f_range[0] / df)
f_max = np.round(f_range[1] / df)

f_vec = np.arange(f_min, f_max+df, 1, dtype=float)
f_vec = f_vec * df
f_vec = f_vec.reshape((f_vec.size,1))

w_vec = 2 * np.pi * f_vec

t_vec = np.arange(0, T, dt, dtype=float)
t_vec = t_vec.reshape((1,t_vec.size))

nf = np.size(f_vec)
nt = np.size(t_vec)

#%% Construct randomized multisine arrays
# Phase randomized between -pi and pi
#np.random.seed(1) # Use to control the random number generator with a seed
ph_mat = -np.pi + (2*np.pi) * (np.random.random_sample((nf,numPhases)))

# Create amplitude vector
Amp = np.ones((np.size(f_vec),1)) / f_vec**c

# Scale amplitude vector to desired RMS
RMS = np.sqrt(np.sum(Amp**2)) / 2**0.5
rms_lim = peak_ampl/2**0.5
gain = (rms_lim / RMS) * atten
Amp = Amp * gain

# Create exp(iwt) and complex amplitude, uses Euler formula
exp_mat = np.exp(1j * np.multiply(w_vec,t_vec))
Amp_f = Amp * np.exp(1j * ph_mat)
signal = np.real(Amp_f.transpose().dot(exp_mat)).transpose()

# Sort in ascending p2p
p2p_diff = signal.max(0) - signal.min(0)
p2p = np.sort(p2p_diff)
ind = np.argsort(p2p_diff)
signal = signal[:,ind]
Amp_f = Amp_f[:,ind]

# Remove any possible signals with p2p larger than 10% the minimum p2p
signal = signal[:,p2p <= p2p[0]*1.1]
Amp_f = Amp_f[:,p2p <= p2p[0]*1.1]

# Generate a list of possible combinations of number of DOFs
temp = np.arange(0,np.size(Amp_f,1))
signal_set = np.array(list(itt.combinations(temp, numDOF)))

#%% Choose the set with the smallest condition to ensure the set is well 
#     conditioned
cur_min = 0
cond_idx = -1
cond_vec = np.zeros((nf,1))

for k in range(np.size(signal_set,0)): #range(signal_set.size):
    for b in range(nf):
        test_signal = Amp_f[b,signal_set[k,:]].reshape((1,numDOF))
        cond_vec[b] = np.linalg.cond(test_signal.dot(test_signal.T))
    if cond_vec.max() < cur_min or cond_idx == -1:
        currMin = cond_vec.max
        cond_idx = k
sig_set_used = np.reshape(signal_set[cond_idx,:],(1,numDOF))

#%% Construct test matrices
for k in range(reps):
    xt = np.tile(signal[:,sig_set_used[k,:]],(reps,1))
    xf = Amp_f[:,sig_set_used[k,:]]
vel = np.gradient(xt,dt, axis=0)
acc = np.gradient(vel,dt, axis=0)

#%% Find frequency spectrum
host_spec = sp.fft.fft(xf)
N = len(host_spec)
n = np.arange(N)
T = N/100
freq = n/T
print(host_spec.shape)
#%% Save multisine in .csv files
if save_test_files:
    test_vals = np.zeros((np.size(xt,0),6))
    test_pos = test_vals
    test_vel = test_vals
    test_acc = test_vals
    ind = np.where(test_DOF!=0)
    for k in range(np.size(xt,1)):
        test_pos[:,ind[0][k]] = xt[:,k]
        test_vel[:,ind[0][k]] = vel[:,k]
        test_acc[:,ind[0][k]] = acc[:,k]
    
    pos_file = save_file_name + " pos" + ".csv"
    vel_file = save_file_name + " vel" + ".csv"
    acc_file = save_file_name + " acc" + ".csv"
    time_file = save_file_name + " time" + ".csv"
    np.savetxt(pos_file,test_pos, delimiter=',')
    np.savetxt(vel_file,test_vel, delimiter=',')
    np.savetxt(acc_file,test_acc, delimiter=',')
    np.savetxt(time_file,t_vec.T, delimiter=',')

#%% Plot multisine signal
units = ["Pos/Rad", "Vel/RadVel", "Acc/RadAcc"]
DOFs = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

fig, axs = plt.subplots(3)
lines = []
print(ind)
for k in range(numDOF):
    axs[k].plot(t_vec.T,xt)
    axs[k].set_ylabel(units[k])
    axs[k].grid(visible=1,which='major',axis='both')
    print(ind[0][k])
    np.append(lines,DOFs[ind[0][k]],axis=0)
lines = DOFs[test_DOF!=0]
print(lines)
plt.legend(lines)
plt.xlabel("Time")

plt.show()
print("Program Complete")