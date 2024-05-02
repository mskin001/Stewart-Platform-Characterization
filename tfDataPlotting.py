import os
from matplotlib import pyplot as plt
import numpy as np

fileName = ["011-Su-p05-p5-p5-p5", "013-Sw-p05-p5-p5-p5", "012-He-p05-p5-p5-p5"]
folder = "Y:\\5700\Water\Marine\LAMP\Projects\\2023 LAMP Characterization\Characterization Data\Transfer Function Data"
dt = 0.01
dofs = ["Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw"]

fig, tfp = plt.subplots(6, 6, sharex=True)
fig.tight_layout()
fig, lagp = plt.subplots(6, 6, sharex=True)
fig.tight_layout()
fig, dof1p = plt.subplots(2, 3, sharex=True)
fig.tight_layout()
fig, dirCmp = plt.subplots(2, 1, sharex=True)
dirCmp[0].set_ylabel("Gain")
dirCmp[1].set_ylabel("Lag [s]")
dirCmp[1].set_xlabel("Freq [Hz]")


for b in range(len(fileName)):
    tfFile = "tf"+fileName[b]+".csv"
    tlFile = "tl"+fileName[b]+".csv"
    fyxFile = "fyx"+fileName[b]+".csv"
    tfPath = os.path.join(folder, tfFile)
    tlPath = os.path.join(folder, tlFile)
    fyxPath = os.path.join(folder, fyxFile)

    tf = np.loadtxt(tfPath, dtype=np.complex64, delimiter=',')
    tl = np.loadtxt(tlPath, dtype=np.complex64, delimiter=',')
    fyx = np.loadtxt(fyxPath, dtype=np.complex64, delimiter=',')

    for k in range(6):
        tfp[k,b].semilogx(fyx, np.abs(tf[:,k]))
        tfp[k,b].grid(visible=1,which="major",axis="both")
        tfp[k,b].vlines(x=0.5, ymin=0, ymax=max(np.abs(tf[:,k])), color="orange")
        tfp[k,0].set_ylabel(dofs[k])
    tfp[0,0].set_xlim((0.005,0.8))
    tfp[0,b].set_title(dofs[b])
    tfp[5,b].set_xlabel("Freq [Hz]")

    for k in range(6):
        lagp[k,b].semilogx(fyx, tl[:,k])
        lagp[k,b].grid(visible=1,which="major",axis="both")
        lagp[k,b].vlines(x=0.5, ymin=0, ymax=max(np.abs(tf[:,k])), color="orange")
        lagp[k,0].set_ylabel(dofs[k])
    lagp[0,b].set_title(dofs[b])
    lagp[5,b].set_xlabel("Freq [Hz]")
    
    dof1p[0,b].semilogx(fyx, np.abs(tf[:,b]))
    dof1p[0,b].grid(visible=1,which="major",axis="both")
    dof1p[1,b].semilogx(fyx,tl[:,b])
    dof1p[1,b].grid(visible=1,which="major",axis="both")
    dof1p[0,b].set_title(dofs[b])
    dof1p[1,b].set_xlabel("Freq [Hz]")

    dirCmp[0].semilogx(fyx, np.abs(tf[:,b]))
    dirCmp[1].semilogx(fyx, np.abs(tl[:,b]))
    dirCmp[0].set_ylabel("Gain")
    dirCmp[1].set_ylabel("Lag [s]")
    dirCmp[1].set_xlabel("Freq [Hz]")
    dirCmp[0].legend(dofs[0:3])

dof1p[0,0].set_ylabel("Gain")
dof1p[1,0].set_ylabel("Lag [s]")

plt.show()
print("Program Complete")