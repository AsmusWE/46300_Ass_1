# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 11:18:09 2024

@author: freja
"""
import numpy as np
from math import pi, exp, sin, cos
import matplotlib.pyplot as plt


#Q2 in Assignment 1

R = 89.17           # m
A = pi * R**2       # m^2
B = 3               # -
P_rated = 10.64     # MW
V_in = 4            # m/s
V_out = 25          # m/s
rho = 1.225         # kg/m^3
Cp_max = 0.46599649          # - INSERT FROM Q1
tip_speed_max = 8   # - INSERT FROM Q1



# Initialize arrays for wind speeds and power output
V_range = np.linspace(V_in, V_out, 500)  # Wind speeds from V_in to V_out
P = np.zeros_like(V_range)               # Power array to store calculated power
omega = np.zeros_like(V_range)           # Omega array to store angular speed
V_rated = None                           # Variable to store wind speed at rated power

# Calculate power for each wind speed and stop when it reaches P_rated
for i, V in enumerate(V_range):
    P[i] = 0.5 * rho * A * V**3 * Cp_max / 1e6  # Convert to MW
    if P[i] >= P_rated and V_rated is None:     # Check if we reached rated power
        V_rated = V                            # Store the wind speed when rated power is first reached
        print(f"Rated power of {P_rated} MW is reached at wind speed: {V_rated:.2f} m/s")
    if P[i] > P_rated:                         # Cap power at rated power
        P[i] = P_rated
    
    # Calculate angular velocity (omega)
    omega[i] = (tip_speed_max * V) / R  # Angular speed in rad/s
    


#Find omega rated
omega_max = (tip_speed_max * V_rated) /R #angular speed, m/s

# Print the rated wind speed and corresponding omega
print(f"The rated power of {P_rated} MW is first reached at a wind speed of {V_rated:.2f} m/s, "
      f"and the maximum angular speed is {omega_max:.2f} rad/s.")

#%%
# Plotting the power curve
plt.figure(figsize=(8, 6))
plt.plot(V_range, P, label='Power Output (MW)', color = 'orange')
plt.axhline(P_rated, color = 'royalblue', linestyle='--', label='Rated Power (10.64 MW)')

plt.title('Wind Turbine Power Curve')
plt.xlabel('Wind Speed (m/s)')
plt.ylabel('Power Output (MW)')
plt.legend()
plt.grid(True)
plt.show()


# Plotting omega as a function of wind speed
plt.figure(figsize=(8, 6))
plt.plot(V_range, omega, label='Angular Velocity (rad/s)', color='green')

plt.title('Angular Velocity vs Wind Speed')
plt.xlabel('Wind Speed (m/s)')
plt.ylabel('Angular Velocity (rad/s)')
plt.grid(True)
plt.legend()
plt.show()

