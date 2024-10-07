import numpy as np
import math

# Defining functions to compute parameters
def compute_flow_angle(a, a_prime, Vo, omega, r):
    denominator = (1 + a_prime) * omega * r
    if np.isclose(denominator, 0):
        return 0  # Avoid division by zero, return a valid angle
    return math.atan((1 - a) * Vo / denominator)  # flow angle phi

def compute_solidity(B, c, r):
    return (B * c) / (2 * math.pi * r)  # solidity sigma

def compute_Vrel(Vo, omega, r):
    return math.sqrt(Vo**2 + (omega * r)**2)

def prandtl_tip_loss_factor(B, R, r, phi):
    sin_phi = np.sin(np.abs(phi))
    exp_term = np.exp(- (B / 2) * (R - r) / (r * sin_phi))
    return (2 / math.pi) * math.acos(exp_term)

def compute_aero_forces(Cl, Cd, phi):
    Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
    Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
    return Cn, Ct

def compute_thrust_coefficient(a, phi, Cn, sigma, F):
    sin_phi = np.sin(phi)
    if np.isclose(F, 0) or np.isclose(sin_phi, 0):
        return 0  # Avoid division by zero
    return ((1 - a)**2 * Cn * sigma) / (F * (sin_phi**2))

def compute_aero_forces(Cl, Cd, phi):
    Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
    Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
    return Cn, Ct

def update_a(a, f, dC_T):
    if a < 1/3:
        a_star = dC_T / (4 * (1 - a))
    else:
        a_star = dC_T / (4 * (1 - 0.25 * a * (5 - 3 * a)))
    return f * a_star + (1 - f) * a

def update_a_prime(a_prime, f, Ct, sigma, phi, F):
    if np.isclose(F, 0) or np.isclose(np.sin(phi) * np.cos(phi), 0):
        return 0  # Avoid division by zero
    a_prime_star = (Ct * sigma / (4 * F * np.sin(phi) * np.cos(phi))) * (1 + a_prime)
    return f * a_prime_star + (1 - f) * a_prime


## ****** INTERPOLATION OF COEFFICIENTS Cl AND Cd ********
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
def force_coeffs_10MW(angle_of_attack,thick,aoa_tab,cl_tab,cd_tab,cm_tab):
    cl_aoa=np.zeros([1,6])
    cd_aoa=np.zeros([1,6])
    cm_aoa=np.zeros([1,6])
    #Interpolate to current angle of attack:
    for i in range(np.size(files)):
        cl_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cl_tab[:,i])
        cd_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cd_tab[:,i])
        cm_aoa[0,i]=np.interp (angle_of_attack,aoa_tab,cm_tab[:,i])
    
#<<<<<<< HEAD
     #Interpolate to current thickness:
    cl=np.interp (thick,thick_prof,cl_aoa[0,:])
    cd=np.interp (thick,thick_prof,cd_aoa[0,:])
    cm=np.interp (thick,thick_prof,cm_aoa[0,:])
    return cl, cd, cm



chords = np.arange(0, 3, 0.5)
theta_p6 =  np.arange(-4, 4, 1)


Vo = 10
r6 = 80.14
R = 89.17 
B = 3
f = 0.1
beta6 = -2.28
tsr6 = 8
omega6 = (tsr6 * Vo) / R  # Calculate omega for the current tip speed ratio
thick6 = 24.10
C_P6_max = 0.0
c_at_max = 0.0
theta_at_max = 0.0
epsilon = 1*10**(-3)

def BEM2(Vo, omega, theta_p, c6):
    a, a_prime = 0.0, 0.0  # Step 1: Re-initialize a and a' for each iteration
    icount = 0  # Initialize iteration counter
    k = 0  # Loop control variable
    
    while k == 0:
        icount = icount + 1
        # Step 2: Compute the flowangle φ from
        phi = compute_flow_angle(a, a_prime, Vo, omega6, r6)
        # Step 3: Compute the local angle of attack by subtracting the twist and global pitch from the flow angle
        alpha = math.degrees(phi) - (beta6 + theta_p)  # Angle of attack
        # Step 4: Interpolate the Cl (α) and Cd(α)
        Cl, Cd, Cm = force_coeffs_10MW(alpha,thick6,aoa_tab,cl_tab,cd_tab,cm_tab)
        # Step 5: Compute Cn and Ct 
        Cn, Ct = compute_aero_forces(Cl, Cd, phi)
        # Step 6: Compute the local thrust coefficient
        F = prandtl_tip_loss_factor(B, R, r6, phi)
        sigma = compute_solidity(B, c6, r6)
        ## Updating a using Glauert correction
        dC_T = compute_thrust_coefficient(a, phi, Cn, sigma, F)
        # Step 7: Update a n and a’ n
        a_new = update_a(a, f, dC_T)
        #a_new = 0.246 * dC_T + 0.0586 * dC_T**2 + 0.0883 * dC_T**3
        a_prime_new = update_a_prime(a_prime, f, Ct, sigma, phi, F)
        # Step 8: If ǀan- an-1ǀ and ǀa’ n- a’ n-1ǀ larger than ε go back and recalculate
        if np.abs(a_new - a) < epsilon  and np.abs(a_prime_new - a_prime) < epsilon:
            k = 1  # Set k to 1 to exit the loop
        else:
            k = 0  # Set k to 0 to cotinue the loop
            
        # Check for maximum iterations
        if icount > 1000:
            #print(f"Not converged for radius {r} with c {c} and beta {beta}.")
            k = 1  # Exit the loop if max iterations reached
        a = a_new
        a_prime = a_prime_new
    

    return a, a_prime, phi, Ct, sigma, F


# Nested for loop to go over each value combination of omega and theta_p
for j, c in enumerate(chords):
    for i, theta in enumerate(theta_p6):
        # Loop through each blade element
        a, a_prime, phi, Ct, sigma, F = BEM2(Vo, omega6, theta, c)  # Assume BEM function returns [pn, pt]
        test = (r6/R) * ((tsr6 * (1-a)**2 * Ct * sigma)/ (np.sin(phi)**2))
        if test > C_P6_max:
            C_P6_max = test
            c_at_max = c
            theta_at_max = theta
            
        
        
print(f"The maximum CP is: {C_P6_max:.3f}, and it occurs at c = {c_at_max:.3f} and pitch = {theta_at_max:.3f} " )