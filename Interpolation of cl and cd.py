# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 10:42:29 2024

@author: freja
"""

import os
import numpy as np


# Define the correct directory where your files are located
file_directory = r'C:\Users\freja\OneDrive\Dokumenter\Freja Skole\1. Semester\43600 mandag+torsdag\Assignment 1\interpolate'

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
    # Create the full file path
    file_path = os.path.join(file_directory, files[i])
    
    # Ensure the file exists
    if os.path.exists(file_path):
        aoa_tab[:], cl_tab[:,i], cd_tab[:,i], cm_tab[:,i] = np.loadtxt(file_path, skiprows=0).T
    else:
        print(f"File {files[i]} not found at {file_path}")


# Thickness of the airfoils considered
# NOTE THAT IN PYTHON THE INTERPOLATION REQUIRES THAT THE VALUES INCREASE IN THE VECTOR!

thick_prof=np.zeros(6)
thick_prof[0]=24.1;
thick_prof[1]=30.1;
thick_prof[2]=36;
thick_prof[3]=48;
thick_prof[4]=60;
thick_prof[5]=100;


def force_coeffs_10MW(angle_of_attack,thick,aoa_tab,cl_tab,cd_tab,cm_tab):
    cl_aoa=np.zeros([1,6])
    cd_aoa=np.zeros([1,6])
    cm_aoa=np.zeros([1,6])
    

    #Interpolate to current angle of attack:
    for i in range(np.size(files)):
        cl_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cl_tab[:,i])
        cd_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cd_tab[:,i])
        cm_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cm_tab[:,i])
    
    #Interpolate to current thickness:
    cl=np.interp (thick,thick_prof,cl_aoa[0,:])
    cd=np.interp (thick,thick_prof,cd_aoa[0,:])
    cm=np.interp (thick,thick_prof,cm_aoa[0,:])


    return cl, cd, cm 



# Lets test it:
angle_of_attack=-10 # in degrees
thick = 27 # in percent !
[clift,cdrag,cmom]=force_coeffs_10MW(angle_of_attack,thick,aoa_tab,cl_tab,cd_tab,cm_tab)

print('cl:',clift)
print('cd:',cdrag)
print('cm:',cmom)