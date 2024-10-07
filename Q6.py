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
    
<<<<<<< HEAD
     #Interpolate to current thickness:
    cl=np.interp (thick,thick_prof,cl_aoa[0,:])
    cd=np.interp (thick,thick_prof,cd_aoa[0,:])
    cm=np.interp (thick,thick_prof,cm_aoa[0,:])
    return cl, cd, cm
=======
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

def V_rel (omega, r, V0):
    """
    Calculate the relative wind speed:
        
    V_0 = relative wind speed 
    omega = omega
    r = radius where we calculate
    """
    V_rel = math.sqrt(V0**2 + (omega * r)**2)
    
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
"""
# Constants and parameters
r = 80.14                   # radius at which we're optimizing (m)
R = 89.17                   # total blade length (m)
B = 3                       # number of blades
V0 = 8                   # free-stream wind speed (m/s)
tip_speed = 8               # tip speed ratio (lambda)
omega = V0 * tip_speed / r  # calculate omega using lambda = 8
beta = -2.28                # beta at the given radius
tc = 24.10                  # thickness cord ratio
f = 0.1                     # Under relaxation needed
rho = 1.225                 # Densit of air (kg/m^3)
A = math.pi * (R**2 - r**2) # Swept rotor area

# Initialize variables for the BEM loop
Ct = 0 
sigma = 0
phi = 0
               
a = 0
a_prime = 0
epsilon = 1e-6  # Set a small threshold for convergence
max_iterations = 1000  # Maximum number of iterations to prevent infinite loops

# Initialize variables
max_Cp = 0  # To store the maximum Cp found
best_c = None           # To store the c value corresponding to the max Cp
best_teta = None       # To store the angle corresponding to the max Cp

# Define the range and step size for c and angle
c_min = 0.1
c_max = 1.7
c_step = 0.01

teta_min = -6# Adjust as per your range of angles
teta_max = -3
teta_step = 0.01

# Loop through c values from c_min to c_max
c = c_min

while c <= c_max:
    
    # Loop through angle values from angle_min to angle_max
    teta = teta_min
    while teta <= teta_max:
                
        #Variable deciding whether or not to use the madsen way or the Glauert way
        madsen = 0    
        
        # Iteration loop
        iteration = 0
        while iteration < max_iterations:
            # Save the current values of a and a_prime
            a_old = a
            a_prime_old = a_prime
            
            # Calculate the flow angle phi
            phi = flow_angle(V0, omega, r, a, a_prime)
            phi_deg = phi*180/math.pi
            
            # Calculate the angle of attack
            alpha_local = angle_of_attack(phi_deg, teta, beta)
    
            # Calculate cl, cd and cm
            C_l, C_d, C_m = InterpolationOfCd.force_coeffs_10MW(alpha_local, tc, aoa_tab, cl_tab, cd_tab, cm_tab)
            
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
                if madsen == 0 :
                    a_temp = dCT / (4 * (1 - 0.25 * (5 - 3 * a_old) * a_old))
                else:
                    a_temp = 0.246 * dCT + 0.0586 * dCT**2 + 0.0883 * dCT**3
            
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

        Cp = (r/R)* ((tip_speed * (1-a)**2 *Ct * sigma) / (math.sin(phi)**2))  # Call your function to calculate Cp based on c and angle
        
        # Check if this Cp is larger than the previous maximum Cp
        if Cp > max_Cp:
            max_Cp = Cp       # Update maximum Cp
            best_c = c        # Update best c corresponding to the maximum Cp
            best_teta = teta  # Update best angle corresponding to the maximum Cp
        
        teta += teta_step  # Increment the angle by the step size
        
    c += c_step  # Increment c by the step size

# Print or return the results
print(f"The maximum Cp is {max_Cp:.4f}, and it occurs at c = {best_c:.4f} and angle = {best_teta:.4f}")
"""

>>>>>>> 277f4ce72fb9f6e6fd2b264ff8aea2dc51632f73


chords = np.arange(0, 3, 0.5)
theta_p6 =  np.arange(-4, 4, 1)


Vo = 10
r6 = 80.14
R = 89.17 
beta6 = -2.28
tsr6 = 8
omega6 = (tsr6 * Vo) / R  # Calculate omega for the current tip speed ratio
thick6 = 24.10
C_P6_max = 0.0
c_at_max = 0.0
theta_at_max = 0.0

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
            print(f"Not converged for radius {r} with c {c} and beta {beta}.")
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