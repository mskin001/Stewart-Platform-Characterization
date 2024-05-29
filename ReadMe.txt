Test Profile File Names:

TPXXX-Su-pXX-pXX-pXX-pXX-HF
^     ^  ^       ^       ^ Ending characters to identify any modifiers (HF - Hex Frame, M1 - Mass 1, etc.)
|     |  |       |-- Amplitdue in the profile with "p" representing the decimal point (p05-1p8 is amplitudes increasing from 0.5m to 1.8m)
|     |  |-- Frequencies in the profile with "p" representing the decimal point (p05-2p0 is frequencies increasing from 0.05Hz to 2.0Hz)
|     |-- DOF Su, Sw, He, Ro, Pi, Ya (Surge, Sway, Heave, Roll, Pitch, Yaw)
|-- Test Profile Number (TP001 = Test Profile 001)

Emulator and Real World file names will follow the same convension but substituting TP for EM or RW for emulator or Real World respectively.

Signal Names:
Data and variables beginning with "h2c" or "host" are signals sent to LAMP from the host. This is what we (the operators) are asking LAMP to
do. Data and variables beginning with "c2h" or "control" or "cont" are signals coming back from LAMP to us. This is the response, what LAMP 
is actually doing.

File Storage:
All data is stored at: Y:\5700\Water\Marine\LAMP\Projects\2023 LAMP Characterization\Characterization Data

Subfolders in this directory organize the Test Profile, Real World, Emulator, etc. files. Data processing files use these directory names to 
load, process, post-process, and save data. The directory names can be changed in lines 8-25 of "process_LAMP_data.py" if needed. 

Files:
ramp_gen.py -> generates a smooth sine wave with smoothing increasing frequency, amplitude, or both. Can produce up to 6 sine waves with
               a prescribed indipendent phase shift. Also produces velocity and amplitude signals using the gradient function. Rotational
               directions (Roll, Pitch, Yaw) require scaling which has not yet been implimented. If selected the output files are saved in 
               the same directory using the given file name.
noise_gen.py -> generates up to 6 indipendent randomly varying multi-sine waves. Amplitude and frequency ranges can be controlled and the 
               output files can be saved in the same directory using the given file name, if selected.
process_LAMP_data.py -> Top level data processing for signal from LAMP's control program "Commander". Reads data from "file_name" and completes
               all the data processing. The resulting transfer functions can be saved with the same file naming convention described above. 
               Data processing functions are imporded from "lampDataFunc.py". Desired output plots are controllable. Each set can produce up to 
               6 unique figures so caution should be used to avoid reaching the plot display limit. Input files MUST have columns in a specific
               order listed below.
lampDataFunc.py -> location for all the data processing functions. New functions should be added here and called in the top level file 
               "process_LAMP_data.py"
tfDataPlotting.py -> post processing file for displaying the transfer function results saved by "process_LAMP_data.py".

LAMP Data File Column Order:
"process_LAMP_data.py" expects the columns in the LAMP data to be in the following order. This can be achieved by reordering the columns in Excel
or creating an additional program to preprocess the raw data files.

Order:
Time, SystemTime, Motion.KinematicsFwd->ComPose0, Motion.KinematicsFwd->ComPoseVel0, Motion.KinematicsFwd->ComPoseAcc0, Motion.KinematicsFwd->ComPose1, Motion.KinematicsFwd->ComPoseVel1, Motion.KinematicsFwd->ComPoseAcc1, Motion.KinematicsFwd->ComPose2, Motion.KinematicsFwd->ComPoseVel2, Motion.KinematicsFwd->ComPoseAcc2, Motion.KinematicsFwd->ComPose3, Motion.KinematicsFwd->ComPoseVel3, Motion.KinematicsFwd->ComPoseAcc3, Motion.KinematicsFwd->ComPose4, Motion.KinematicsFwd->ComPoseVel4, Motion.KinematicsFwd->ComPoseAcc4, Motion.KinematicsFwd->ComPose5, Motion.KinematicsFwd->ComPoseVel5, Motion.KinematicsFwd->ComPoseAcc5, Motion.KinematicsBwd->Pose0, Motion.KinematicsBwd->Pose1, Motion.KinematicsBwd->Pose2, Motion.KinematicsBwd->Pose3, Motion.KinematicsBwd->Pose4, Motion.KinematicsBwd->Pose5, Motion.KinematicsBwd->PoseVel0, Motion.KinematicsBwd->PoseVel1, Motion.KinematicsBwd->PoseVel2, Motion.KinematicsBwd->PoseVel3, Motion.KinematicsBwd->PoseVel4, Motion.KinematicsBwd->PoseVel5, 
Motion.KinematicsBwd->PoseAcc1, Motion.KinematicsBwd->PoseAcc2, Motion.KinematicsBwd->PoseAcc3, Motion.KinematicsBwd->PoseAcc4, Motion.KinematicsBwd->PoseAcc5

Coresponding Names:
Time, SystemTime, H2C Surge Pos, H2C Surve Vel, H2C Surge Acc, H2C Sway Pos, H2C Sway Vel, H2C Sway Acc, H2C Heave Pos, H2C Heave Vel, H2C Heave Acc, H2C Roll Pos, H2C Roll Vel, H2C Roll Acc, H2C Pitch Pos, H2C Pitch Vel, H2C Pitch Acc, H2C Yaw Pos, H2C Yaw Vel, H2C Yaw Acc, C2H Surge Pos, C2H Sway Pos, C2H Heave Pos, C2H Roll Pos, C2H Pitch Pos, C2H Yaw Pos, C2H Surge Vel, C2H Sway Vel, C2H Heave Vel, C2H Roll Vel, C2H Pitch Vel, C2H Yaw Vel, C2H Surge Acc, C2H Sway Acc, C2H Heave Acc, C2H Roll Acc, C2H Pitch Acc, C2H Yaw Acc

(Look I know the order is totally asinine. One of the early output files from the emulator had this column order and I just never changed it. It's a vestigial relic from the creation of this program. Don't H8 me.)