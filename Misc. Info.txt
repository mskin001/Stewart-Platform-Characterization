Test Profile File Names:

TP000-000000-FFF-A-PVA
^     ^      ^   ^-- Amplitdue in the profile in [m]
|     |      |-- Frequencies in the profile, last digit is the largest frequency in Hz (052 = 0.05Hz - 2 Hz)
|     |-- DOF in the test profile with 1 included and 0 not included (Surge, Sway, Heave, Roll, Pitch, Yaw) (100010 = Surge, Pitch)
|-- Test Profile Number (TP001 = Test Profile 001)

Emulator and Real World file names will follow the same convension but
substituting TP for EM or RW for emulator or Real World respectively.

Signal Names:
Data and variables beginning with "h2c" or "host" are signals sent to
LAMP from the host. This is what we (the operators) are asking LAMP to
do. Data and variables beginning with "c2h" or "control" or "cont" 
are signals coming back from LAMP to us. This is the response, what
LAMP is actually doing.