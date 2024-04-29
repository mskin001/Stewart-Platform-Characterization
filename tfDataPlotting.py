import os
from matplotlib import pyplot as plt
import numpy as np

fileName = ["011-Su-p05-p5-p5-p5", "013-Sw-p05-p5-p5-p5", "012-He-p05-p5-p5-p5"]
folder = "Y:\\5700\Water\Marine\LAMP\Projects\\2023 LAMP Characterization\Characterization Data\Transfer Function Data"
dt = 0.01

fig, tfp = plt.subplots(6, 6, sharex=True)
plt.tight_layout()
fig, lagp = plt.subplots(6, 6, sharex=True)
plt.tight_layout()

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
        tfp[k,b].semilogx(fyx, (np.abs(tf[:,k])))
        tfp[k,b].set_ylabel("Gain")
        tfp[k,b].grid(visible=1,which="major",axis="both")
    tfp[0,0].set_xlim((0.005,0.8))

    for k in range(6):
        lagp[k,b].semilogx(fyx, tl[:,k])
        lagp[k,b].grid(visible=1,which="major",axis="both")
        lagp[k,b].set_ylabel("Lag [sec]")
    

plt.show()
print("Program Complete")