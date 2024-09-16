# -*- coding: utf-8 -*-
"""
Spyder Editor

Created on Mon Sep 9 09:38:31 2024

@author: Freja
"""
import numpy
import math


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

def thrust_coefficient (C_n, sigma, F, phi):
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

def p_n (rho, Vrel, c, Cn):
    """
    Calculate the Normal force per unit length
    
    p_n = normal force per meter
    rho = density of the air
    c = cord length
    C_n = Normal force coefficient
    """
    p_n = 1/2 * rho * Vrel**2 * c * Cn
    
    return p_n

def p_t (rho, Vrel, c, Ct):
    """
    Calculate the tangential force per unit length
    
    p_t = tangential force per meter
    rho = density of the air
    c = cord length
    C_t = tangential force coefficient
    """
    p_t = 1/2 * rho * Vrel**2 * c * Ct
    
    return p_t


#Now the while loop calculating the error in a^n and a^n-1 and the same for a_prime
# Initialize variables
a = 0
a_prime = 0
epsilon = 1e-10  # Set a small threshold for convergence
max_iterations = 1000  # Maximum number of iterations to prevent infinite loops

# Define known parameters (you need to specify these values)
r = 24.5                    # Radius                   (m)
c = 1.5                     # cord length              (m)

R = 31                      # Full blade ratio         (m)
B = 3                       # Number of blades
rho = 1.225                 # Densit of air            (kg/m^3)
V0 = 8.0                    # Example wind speed       (m/s)
omega = 2.61                    # Example angular velocity (rad/s)
teta = -3.0                 # Global pitch angle       (degree)
beta = 2.0                  # Twist angle              (radians)
C_l = 0.5                   # Lift coefficient 
C_d = 0.01                  # Drag coefficient 
f = 0.1                     # Under relaxation needed


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
    dCT = thrust_coefficient(Cn, sigma, F, phi)
    
    
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
        print(f"Converged after {iteration} iterations.")
        break
    
    iteration += 1
else:
    print("Maximum iterations reached without convergence.")

# Final values of a and a_prime, rounded to 3 decimals
print(f"a = {a:.3f}, a_prime = {a_prime:.3f}")



# Calculate the relative wind speed
Vrel = V_rel(omega, r, V0)

# Calculate the normal and tangential force per meter
pn = p_n(rho, Vrel, c, Cn)
pt = p_t(rho, Vrel, c, Ct) 

print(f"pn = {pn:.3f}, pt = {pt:.3f}")
    

