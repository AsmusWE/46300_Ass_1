# -*- coding: utf-8 -*-
"""
Spyder Editor

Created on Mon Sep 9 09:38:31 2024

@author: Freja
"""
import numpy as np
import math
import matplotlib.pyplot as plt


def flow_angle(V0, omega, r, a, a_prime):
    """
    Calculate the flow angle:
        
    phi = flow angle 
    V_0 = velocity of the wind
    w = omega
    r = radius
    """
    phi = math.atan(((1 - a) * V0) / ((1 + a_prime) * omega * r))
    
    return phi


def angle_of_attack (phi, teta, beta):
    """
    Calculate the angle of attack:
        
    alpha local = local angle of attack
    phi = flow angle 
    Teta = global pitch
    beta = twist angle
    """
    alpha_local = phi - (teta + beta)
    
    return alpha_local


def C_n (C_l, C_d, phi):
    """
    Calculate C_n:
        
    C_n = normal force coefficient
    C_l = lift coefficient
    C_d = drag coefficient
    phi = flow angle 
    """
    C_n = (C_l * math.cos(phi)) + (C_d * math.sin(phi))
    
    return C_n


def C_t (C_l, C_d, phi):
    """
    Calculate C_t:
        
    C_t = tangential force coefficient
    C_l = lift coefficient
    C_d = drag coefficient
    phi = flow angle 
    """
    C_t = (C_l * math.sin(phi)) - (C_d * math.cos(phi))
    
    return C_t

def thrust_coefficient (C_n, sigma, F, phi,a):
    """
    Calculate C_T:
        
    dC_T = Thrust force coefficient
    C_n = normal force coefficient
    sigma = solidity ratio
    phi = flow angle
    """
    dC_T = ((1 - a)**2 * C_n * sigma) / (F * math.sin(phi)**2)
    
    return(dC_T)

def tip_loss_factor (B, R, r, phi):
    """
    Calculate F:
        
    F = Prandtl's tip loss factor
    B = number of blades
    R = full blade ratio
    r = current radial position
    phi = flow angle 
    """
    F = (2 / math.pi) * math.acos(math.exp(-(B / 2) * ((R - r) / (r * math.sin(abs(phi))))))
    
    return F

def update_a (F, phi, sigma, C_n):
    """
    Calculate a:
        
    a = axial induction factor
    F = Prandtl's tip loss factor
    phi = flow angle
    sigma = solidity ratio
    C_n = normal force coefficient
    """
    a = 1 / (((4 * F * math.sin(phi)**2) / (sigma * C_n)) + 1)
    
    return a


def update_a_prime (F, phi, sigma, C_t):
    """
    Calculate a':
        
    a' = tangential induction factor
    F = Prandtl's tip loss factor
    phi = flow angle
    sigma = solidity ratio
    C_t = tangential force coefficient
    """
    a_prime = 1 / (((4 * F * math.sin(phi) * math.cos(phi)) / (sigma * C_t))-1)
    
    return a_prime

def solidity (c, B, r):
    """
    Calculate solidity(sigma):
        
    c = cord length of the blade
    B = number of blades
    r = radius where we calculate
    """
    sigma = (c * B) / (2 * math.pi * r)
    
    return sigma

def V_rel (omega, r, V_0):
    """
    Calculate the relative wind speed:
        
    V_0 = relative wind speed 
    omega = omega
    r = radius where we calculate
    """
    V_rel = math.sqrt(V_0**2 + (omega * r)**2)
    
    return V_rel

def interpol_c(interpol_r):
    r = [2.80, 11.00, 16.87, 22.96, 32.31, 41.57, 50.41, 58.53, 65.75, 71.97, 77.19, 78.71, 80.14, 82.71, 84.93, 86.83, 88.45, 89.17]
    c = [5.38, 5.45, 5.87, 6.18, 6.02, 5.42, 4.70, 4.00, 3.40, 2.91, 2.54, 2.43, 2.33, 2.13, 1.90, 1.63, 1.18, 0.60]

    if interpol_r <= r[0]:
        return c[0]
    elif interpol_r >= r[-1]:
        return c[-1]
    else:
        for i in range(len(r) - 1):
            if r[i] <= interpol_r <= r[i + 1]:
                # Linear interpolation
                return c[i] + (c[i + 1] - c[i]) * (interpol_r - r[i]) / (r[i + 1] - r[i])
            
def interpol_beta(interpol_r):
    r = [2.80, 11.00, 16.87, 22.96, 32.31, 41.57, 50.41, 58.53, 65.75, 71.97, 77.19, 78.71, 80.14, 82.71, 84.93, 86.83, 88.45, 89.17]
    beta = [14.50, 14.43, 12.55, 8.89, 6.38, 4.67, 2.89, 1.21, -0.13, -1.11, -1.86, -2.08, -2.28, -2.64, -2.95, -3.18, -3.36, -3.43]
    
    if interpol_r <= r[0]:
        return beta[0]
    elif interpol_r >= r[-1]:
        return beta[-1]
    else:
        for i in range(len(r) - 1):
            if r[i] <= interpol_r <= r[i + 1]:
                # Linear interpolation
                return beta[i] + (beta[i + 1] - beta[i]) * (interpol_r - r[i]) / (r[i + 1] - r[i])

 

def calc_loads(teta,tip_s_ratio,r,V0):
    """
    Calculate the power coefficient:
        
    lambda = tip speed ratio
    theta = pitch angle
    data = data from the turbine
    """
    # Initialize variables
    a = 0
    a_prime = 0
    epsilon = 1e-10  # Set a small threshold for convergence
    max_iterations = 1000  # Maximum number of iterations to prevent infinite loops

    # Define known parameters
    #r = 24.5                    # Radius                   (m)
    #c = 0.5                     # cord length              (m)
    

    R = 89.17                   # Full blade ratio         (m)
    B = 3                       # Number of blades
    rho = 1.225                 # Densit of air            (kg/m^3)
    #V0 = 10  
    omega = tip_s_ratio*V0/R                 # Example wind speed       (m/s)
    #omega = 2.61                    # Example angular velocity (rad/s)
    #teta = -3.0                 # Global pitch angle       (degree)
    #beta = 2.0                  # Twist angle              (radians)
    C_l = 0.5                   # Lift coefficient 
    C_d = 0.01                  # Drag coefficient 
    f = 0.1                     # Under relaxation needed

    beta = interpol_beta(r)
    c = interpol_c(r)
    # Iteration loop
    iteration = 0
    while iteration < max_iterations:
        # Save the current values of a and a_prime
        a_old = a
        a_prime_old = a_prime
        
        # Calculate the flow angle phi
        phi = flow_angle(V0, omega, r, a, a_prime)
        
        # Calculate the angle of attack
        alpha_local = angle_of_attack(phi, teta, beta)
        
        # Calculate the solidity
        sigma = solidity(c, B, r)
        
        # Calculate C_n and C_t
        Cn = C_n(C_l, C_d, phi)
        Ct = C_t(C_l, C_d, phi)
        
        # Calculate Prandtl's tip loss factor
        F = tip_loss_factor(B, R, r, phi)
        
        # Calculate the thrust coefficient 
        dCT = thrust_coefficient(Cn, sigma, F, phi, a)
        
        
        # Apply Glauert correction for a
        if a_old < 1/3:
            a_temp = dCT / (4 * (1 - a_old))
        else:
            a_temp = dCT / (4 * (1 - 0.25 * (5 - 3 * a_old) * a_old))
        
        # Apply under-relaxation
        a = f * a_temp + (1 - f) * a_old
        
        # Apply correction for a_prime
        a_prime_temp = (Ct * sigma / (4 * F * math.sin(phi) * math.cos(phi))) * (1 + a_prime_old)
        a_prime = f * a_prime_temp + (1 - f) * a_prime_old
        
        # Check convergence
        error_a = abs(a - a_old)
        error_a_prime = abs(a_prime - a_prime_old)
        
        # If both errors are below the threshold, stop the loop
        if error_a < epsilon and error_a_prime < epsilon:
            #print(f"Converged after {iteration} iterations.")
            break
        
        iteration += 1
    else:
        print("Maximum iterations reached without convergence.")
    
    # Compute V_rel
    V_rel = math.sqrt((omega*r+a_prime*omega*r)**2 + (V0-a*V0)**2)

    # Compute loads
    P_n = 1/2*rho*V_rel**2*c*Cn
    P_t = 1/2*rho*V_rel**2*c*Ct

    return P_n, P_t

if __name__ == "__main__":
    V0 = 10
    R = 89.17 # This should probably be done in some nice way where it isn't defined in both function and main
    B = 3 # This too
    th_step = 1
    th_max = 3+th_step
    ti_step = 1
    ti_max = 10+ti_step
    r_step = 1
    theta = np.arange(-4, th_max, th_step)
    tip_s_ratio = np.arange(5, ti_max, ti_step)
    r = np.arange(3, 90, r_step)

    Cp_array = np.zeros((len(theta), len(tip_s_ratio)))
    Ct_array = np.zeros((len(theta), len(tip_s_ratio)))

    Cp_coords = np.zeros((len(theta), len(tip_s_ratio), 2))
    Ct_coords = np.zeros((len(theta), len(tip_s_ratio), 2))

    for th_count, th in enumerate(theta):
        for ti_count, ti in enumerate(tip_s_ratio):
            omega = ti*V0/R
            P_sum = 0
            T_sum = 0
            for ra in r:
                # Calculating loads based on the given parameters
                P_n, P_t = calc_loads(th, ti, ra, V0)
                P_sum += omega*B*P_t*ra*r_step
                T_sum += B * P_n * r_step
            Cp_array[th_count, ti_count] = P_sum
            Ct_array[th_count, ti_count] = T_sum
            print(f"Location: {th_count, ti_count} out of {len(theta), len(tip_s_ratio)}")

    fig,ax = plt.subplots()
    CS = ax.contour(tip_s_ratio, theta, Cp_array)
            




                
    
    
    