from Q1_BEM import *
import matplotlib.pyplot as plt

if __name__ == "__main__":

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


    V0 = 10
    R = 89.17 # This should probably be done in some nice way where it isn't defined in both function and main
    B = 3 # This too
    r_step = 3
    r_max = 89.17+r_step-1
    r = np.arange(3, 90, r_step)

    #Variable deciding whether or not to use the madsen way or the Glauert way
    madsen = 1
    
    # Data for going from loads to cp/cn
    rho = 1.225
    A = R**2 * math.pi

    ti = 8 # Tip speed ratio
    P_max = 10.64*10**6
    #P_max = P_max*10
    omega_max = 1.04 
    V0_list = np.arange(10,26,0.1)
    #V0_list = [11.64,11.64]
    Theta_list = np.arange(0,len(V0_list),dtype=float)
    th = 0
    th_step = 0.1

    madsen = 0

    loops = 0
    max_loops = 100
    for count,V0 in enumerate(V0_list):
        loops = 0
        while loops < max_loops:
            omega = ti*V0/R
            P_sum = 0
            T_sum = 0
            for ra in r:
                # Calculating loads based on the given parameters
                P_n, P_t, a = calc_loads(th, ti, ra, V0, cl_tab, cd_tab, cm_tab, aoa_tab, madsen)
                P_sum += omega*B*P_t*ra*r_step
                T_sum += B * P_n * r_step
            P = P_sum
            
            if P>P_max:
                th += th_step
                loops += 1
                print(f'Loop: {loops}, th: {th}, p: {P}')
            else:
                Theta_list[count-1] = float(th)
                print(f'Count: {count}, out of {len(V0_list)}')
                break
    
    plt.figure()
    plt.plot(V0_list, Theta_list, label='Theta vs V0')
    plt.xlabel('V0 (m/s)')
    plt.ylabel('Theta (degrees)')
    plt.title('Theta vs V0')
    plt.legend()
    plt.grid(True)
    plt.show()

