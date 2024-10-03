import pandas as pd
from Q1_BEM import calc_loads
import numpy as np

R = 89.17
tip_s_ratio = 8
V0 = [5,9,11,20]
B = 3
# Loading data from ashes
Thrust5_ash = pd.read_csv('Ashes Data/5Thrust.txt', delimiter=';', decimal=',', thousands='.', header=None).T
Torque5_ash = pd.read_csv('Ashes Data/5Torque.txt', delimiter=';', decimal=',', thousands='.', header=None).T
Thrust9_ash = pd.read_csv('Ashes Data/9Thrust.txt', delimiter=';', decimal=',', thousands='.', header=None).T
Torque9_ash = pd.read_csv('Ashes Data/9Torque.txt', delimiter=';', decimal=',', thousands='.', header=None).T
Thrust11_ash = pd.read_csv('Ashes Data/11Thrust.txt', delimiter=';', decimal=',', thousands='.', header=None).T
Torque11_ash = pd.read_csv('Ashes Data/11Torque.txt', delimiter=';', decimal=',', thousands='.', header=None).T
Thrust20_ash = pd.read_csv('Ashes Data/20Thrust.txt', delimiter=';', decimal=',', thousands='.', header=None).T
Torque20_ash = pd.read_csv('Ashes Data/20Torque.txt', delimiter=';', decimal=',', thousands='.', header=None).T
PowerThrustDF = pd.read_csv('Ashes Data/PowerThrust.txt', delimiter=';', decimal=',', thousands='.', header=None).T

# Calculating loads 
# START OF TERRIBLE IMPORTED INTERPOLATE PART
current_dir = os.getcwd()  # gets the path of the folder where the code is run

# List of files
files = ['FFA-W3-2411.txt', 'FFA-W3-301.txt', 'FFA-W3-360.txt', 'FFA-W3-480.txt', 'FFA-W3-600.txt', 'cylinder.txt']

# Initialize the arrays (assuming you already know the sizes of aoa_tab, cl_tab, cd_tab, cm_tab)
#Initializing tables    
cl_tab=np.zeros([105,6])
cd_tab=np.zeros([105,6])
cm_tab=np.zeros([105,6])
aoa_tab=np.zeros([105,])

# Read in the tables once at startup of simulation
for i in range(np.size(files)):
    file_path = os.path.join(current_dir, files[i])  # Create the full file path
    aoa_tab[:], cl_tab[:,i], cd_tab[:,i], cm_tab[:,i] = np.loadtxt(file_path, skiprows=0).T


# Thickness of the airfoils considered
# NOTE THAT IN PYTHON THE INTERPOLATION REQUIRES THAT THE VALUES INCREASE IN THE VECTOR!

thick_prof=np.zeros(6)
thick_prof[0]=24.1;
thick_prof[1]=30.1;
thick_prof[2]=36;
thick_prof[3]=48;
thick_prof[4]=60;
thick_prof[5]=100;

# END OF TERRIBLE IMPORTED INTERPOLATE PART

Theta_list = [0,0,0,11.1]
tip_s_ratio = 8
V0_list = [5,9,11,20]
r_values = np.linspace(0, R, len(Thrust5_ash))
Thrust_calc = np.zeros((len(r_values),4))
Torque_calc = np.zeros((len(r_values),4))
r_step = r_values[1] - r_values[0]
P_sum = np.zeros(4)
T_sum = np.zeros(4)

# Calculating load curves for each wind speed
for i in range(4):
    omega = tip_s_ratio*V0[i]/R
    for count, r in enumerate(r_values):
        Thrust_calc[count,i], Torque_calc[count,i], a = calc_loads(Theta_list[i], tip_s_ratio, r, V0_list[i], cl_tab, cd_tab,cm_tab, aoa_tab, 1)
        
        # Checking for nan values
        if Thrust_calc[count,i] == Thrust_calc[count,i]:
            T_sum[i] += B * Thrust_calc[count,i] * r_step
            P_sum[i] += omega*B*Torque_calc[count,i]*r*r_step
Thrust_calc = np.nan_to_num(Thrust_calc)
Torque_calc = np.nan_to_num(Torque_calc)

import matplotlib.pyplot as plt

# Plotting Thrust Comparison in a 2x2 subplot
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Thrust at 5 m/s
axs[0, 0].plot(r_values, Thrust_calc[:, 0], label='Calculated Thrust at 5 m/s')
axs[0, 0].plot(r_values, Thrust5_ash.iloc[:, 0], label='Ashes Thrust at 5 m/s')
axs[0, 0].set_xlabel('Radius (r)')
axs[0, 0].set_ylabel('Thrust')
axs[0, 0].legend()
axs[0, 0].set_title('Thrust Comparison at 5 m/s')

# Thrust at 9 m/s
axs[0, 1].plot(r_values, Thrust_calc[:, 1], label='Calculated Thrust at 9 m/s')
axs[0, 1].plot(r_values, Thrust9_ash.iloc[:, 0], label='Ashes Thrust at 9 m/s')
axs[0, 1].set_xlabel('Radius (r)')
axs[0, 1].set_ylabel('Thrust')
axs[0, 1].legend()
axs[0, 1].set_title('Thrust Comparison at 9 m/s')

# Thrust at 11 m/s
axs[1, 0].plot(r_values, Thrust_calc[:, 2], label='Calculated Thrust at 11 m/s')
axs[1, 0].plot(r_values, Thrust11_ash.iloc[:, 0], label='Ashes Thrust at 11 m/s')
axs[1, 0].set_xlabel('Radius (r)')
axs[1, 0].set_ylabel('Thrust')
axs[1, 0].legend()
axs[1, 0].set_title('Thrust Comparison at 11 m/s')

# Thrust at 20 m/s
axs[1, 1].plot(r_values, Thrust_calc[:, 3], label='Calculated Thrust at 20 m/s')
axs[1, 1].plot(r_values, Thrust20_ash.iloc[:, 0], label='Ashes Thrust at 20 m/s')
axs[1, 1].set_xlabel('Radius (r)')
axs[1, 1].set_ylabel('Thrust')
axs[1, 1].legend()
axs[1, 1].set_title('Thrust Comparison at 20 m/s')

# Adjust layout
plt.tight_layout()
plt.show()

# Plotting Torque Comparison in a 2x2 subplot
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Torque at 5 m/s
axs[0, 0].plot(r_values, Torque_calc[:, 0], label='Calculated Torque at 5 m/s')
axs[0, 0].plot(r_values, Torque5_ash.iloc[:, 0], label='Ashes Torque at 5 m/s')
axs[0, 0].set_xlabel('Radius (r)')
axs[0, 0].set_ylabel('Torque')
axs[0, 0].legend()
axs[0, 0].set_title('Torque Comparison at 5 m/s')

# Torque at 9 m/s
axs[0, 1].plot(r_values, Torque_calc[:, 1], label='Calculated Torque at 9 m/s')
axs[0, 1].plot(r_values, Torque9_ash.iloc[:, 0], label='Ashes Torque at 9 m/s')
axs[0, 1].set_xlabel('Radius (r)')
axs[0, 1].set_ylabel('Torque')
axs[0, 1].legend()
axs[0, 1].set_title('Torque Comparison at 9 m/s')

# Torque at 11 m/s
axs[1, 0].plot(r_values, Torque_calc[:, 2], label='Calculated Torque at 11 m/s')
axs[1, 0].plot(r_values, Torque11_ash.iloc[:, 0], label='Ashes Torque at 11 m/s')
axs[1, 0].set_xlabel('Radius (r)')
axs[1, 0].set_ylabel('Torque')
axs[1, 0].legend()
axs[1, 0].set_title('Torque Comparison at 11 m/s')

# Torque at 20 m/s
axs[1, 1].plot(r_values, Torque_calc[:, 3], label='Calculated Torque at 20 m/s')
axs[1, 1].plot(r_values, Torque20_ash.iloc[:, 0], label='Ashes Torque at 20 m/s')
axs[1, 1].set_xlabel('Radius (r)')
axs[1, 1].set_ylabel('Torque')
axs[1, 1].legend()
axs[1, 1].set_title('Torque Comparison at 20 m/s')

# Adjust layout
plt.tight_layout()
plt.show()

# Plotting Thrust vs V0
plt.figure()
plt.plot(V0, T_sum/1000, label='Calculated Thrust')
plt.plot(V0, PowerThrustDF.iloc[1, :], label='Ashes Thrust')
plt.xlabel('Wind Speed (V0)')
plt.ylabel('Thrust')
plt.legend()
plt.title('Thrust vs Wind Speed')
plt.show()

# Plotting Torque vs V0
plt.figure()
plt.plot(V0, P_sum/1000, label='Calculated Power')
plt.plot(V0, PowerThrustDF.iloc[0, :], label='Ashes Power')
plt.xlabel('Wind Speed (V0)')
plt.ylabel('Torque')
plt.legend()
plt.title('Torque vs Wind Speed')
plt.show()

