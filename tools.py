import numpy as np   # 2.0.2
import pandas as pd  # 2.2.3
import os

# Define Functions

def Read_Data( CASE ):
    # Check File path
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        base_dir = os.getcwd()

    data_dir = os.path.join(base_dir, "CASE_" + CASE + "_INPUT")

    required_files = [
        "CASE_" + CASE + "_BUS DATA.csv",
        "CASE_" + CASE + "_BRANCH DATA.csv",
        "CASE_" + CASE + "_GEN DATA.csv",
    ]

    for f_name in required_files:
        file_path = os.path.join(data_dir, f_name)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Missing Input File: {f_name}")

    # Read Files
    Bus_Data = (
        pd.read_csv(os.path.join(data_dir, "CASE_" + CASE + "_BUS DATA.csv"))
        .to_numpy()
    )
    Bus_Data = Bus_Data[Bus_Data[:, 0].argsort()]  # Make sure bus numbers ascendingly ordered

    # Check if Bus number is sequentially alocated or not
    bus_ids = Bus_Data[:, 0].astype(int)
    expected_ids = np.arange(1, len(bus_ids) + 1)
    if not np.array_equal(bus_ids, expected_ids):
        raise ValueError("Bus Numbers are not sequential (1 to N). Re-indexing is required.")

    Branch_Data = pd.read_csv(
        os.path.join(data_dir, "CASE_" + CASE + "_BRANCH DATA.csv")
    ).to_numpy()
    Gen_Data = pd.read_csv(
        os.path.join(data_dir, "CASE_" + CASE + "_GEN DATA.csv")
    ).to_numpy()

    # # If you need sanity check,
    # print("Bus_Data shape:", Bus_Data.shape)
    # print(Bus_Data)
    # print("Branch_Data shape:", Branch_Data.shape)
    # print(Branch_Data)
    # print("Gen_Data shape:", Gen_Data.shape)
    # print(Gen_Data)

    return Bus_Data, Branch_Data, Gen_Data


def admittance_matrix(N_bus, Bus_Data, Branch_Data, S_base):
    Y = np.zeros( (N_bus, N_bus), dtype=complex )

    N_branch = Branch_Data.shape[0]
    for i in range(N_branch):
        from_bus = int( Branch_Data[i, 0] )
        to_bus   = int( Branch_Data[i, 1] )

        if max(from_bus, to_bus) > N_bus:
            raise ValueError(f"Branch data contains bus no. {max(from_bus, to_bus)} which exceeds N_bus={N_bus}.")

        z_series     = ( Branch_Data[i, 2] ) + ( Branch_Data[i, 3] ) * 1j   # R + jX
        y_shunt_half = ( Branch_Data[i, 4] / 2 ) * 1j                       # j(B/2)

        Y[from_bus-1, from_bus-1] += 1 / z_series + y_shunt_half
        Y[to_bus-1, to_bus-1]     += 1 / z_series + y_shunt_half

        Y[from_bus-1, to_bus-1]   += -1 / z_series
        Y[to_bus-1, from_bus-1]   += -1 / z_series

    for i in range(N_bus):
        Gs_pu = Bus_Data[i, 4] / S_base
        Bs_pu = Bus_Data[i, 5] / S_base

        Y[i, i] -= (Gs_pu) + (Bs_pu) * 1j        # Account for Self Shunt Power

    return Y



def Approx_starting_point(Y_bus, Bus_Data, Gen_Data, S_base, N_bus, N_pq, N_pv, map_tuple, dim, iterations=5):
    """
    Use Fast Decoupled Power Flow (FDPF) logic to generate a better initial guess.
    Assumptions:
    1. Conductance G is negligible (or R << X).
    2. cos(theta) ~ 1, sin(theta) ~ 0.
    """
    
    # 1. Initialize V, theta (Same logic as Flat Start)
    V = np.zeros( (N_bus, 2) )
    gen_v_map = dict(zip(Gen_Data[:, 0], Gen_Data[:, 3]))

    for i in range(N_bus):
        bus_num  = int(Bus_Data[i, 0])
        bus_type = Bus_Data[i, 1]

        if bus_type == 3:              # Slack
            if bus_num in gen_v_map:
                V[i, 0] = gen_v_map[bus_num] 
            else:
                V[i, 0] = 1.0
            V[i, 1] = np.deg2rad(Bus_Data[i, 7]) 

        elif bus_type == 2:            # PV
            if bus_num in gen_v_map:
                V[i, 0] = gen_v_map[bus_num] 
            else:
                V[i, 0] = 1.0 # Fallback
            V[i, 1] = 0.0

        elif bus_type == 1:            # PQ
            V[i, 0] = 1.0
            V[i, 1] = 0.0

    # 2. Construct Constant B Matrices (Approximation using Imaginary part of Y_bus)
    # B' matrix for P-theta (size: N_pq+N_pv x N_pq+N_pv)
    # Uses indices of all non-slack buses
    idx_p_theta = np.array(map_tuple[:(N_pq+N_pv)], dtype=int) - 1
    B_prime = Y_bus.imag[np.ix_(idx_p_theta, idx_p_theta)]

    # B'' matrix for Q-V (size: N_pq x N_pq)
    # Uses indices of PQ buses only
    # Note: In map_tuple, the last N_pq elements are the PQ bus numbers for Voltage variables
    idx_q_v = np.array(map_tuple[(N_pq+N_pv):], dtype=int) - 1
    B_double_prime = Y_bus.imag[np.ix_(idx_q_v, idx_q_v)]

    # 3. FDPF Iteration Loop
    print(f" > Calculating Approximation (FDPF {iterations} iters)...")
    
    for k in range(iterations):
        # Calculate Mismatch (P and Q)
        f = Mismatch(map_tuple, Y_bus, V, Bus_Data, Gen_Data, N_bus, dim, N_pq, N_pv, S_base)
        
        # --- P-Theta Step ---
        # dP = - B' * dTheta  =>  dTheta = - inv(B') * (dP / V)
        dP = f[:(N_pq+N_pv)]
        V_theta_buses = V[idx_p_theta, 0]
        
        # Normalize Mismatch by Voltage Magnitude (Standard FDPF practice)
        rhs_p = dP / V_theta_buses
        
        # Solve for angle correction
        dTheta = - np.linalg.solve(B_prime, rhs_p)
        
        # Update Angles
        V[idx_p_theta, 1] += dTheta

        # --- Q-V Step ---
        # Re-calculate Mismatch with updated angles (Optional but more accurate)
        f_mid = Mismatch(map_tuple, Y_bus, V, Bus_Data, Gen_Data, N_bus, dim, N_pq, N_pv, S_base)
        
        dQ = f_mid[(N_pq+N_pv):]
        V_pq_buses = V[idx_q_v, 0]
        
        # Normalize Mismatch by Voltage Magnitude
        rhs_q = dQ / V_pq_buses
        
        # Solve for magnitude correction
        dV = - np.linalg.solve(B_double_prime, rhs_q)
        
        # Update Magnitudes
        V[idx_q_v, 0] += dV

    # 4. Construct Result Vector x0
    x0 = np.zeros(dim)
    # Angles (PQ + PV)
    x0[:(N_pq+N_pv)] = V[idx_p_theta, 1]
    # Magnitudes (PQ)
    x0[(N_pq+N_pv):] = V[idx_q_v, 0]

    # print(f"Initial Guess: {x0}")

    return x0



def Jacobian( map_tuple, Y_bus, V, N_bus, dim, N_pq, N_pv ):

    if len(map_tuple) != dim:
        raise ValueError(f"Dimension Mismatch !")

    J = np.zeros( (dim, dim) )
    map_idx = np.array( map_tuple, dtype=int ) - 1

    # 1) Real, differentiated by voltage angle (J11)
    for i, bus_idx_i in enumerate(map_idx[:(N_pq+N_pv)]):
        for j, bus_idx_j in enumerate(map_idx[:(N_pq+N_pv)]):
            if bus_idx_i == bus_idx_j:
                for k in range(N_bus):
                    if k == bus_idx_i: continue
                    G_ik     = Y_bus[bus_idx_i, k].real
                    B_ik     = Y_bus[bus_idx_i, k].imag
                    V_i      = V[bus_idx_i, 0]
                    V_k      = V[k, 0]
                    theta_ik = V[bus_idx_i, 1] - V[k, 1]
                    J[i, j] += V_i * V_k * ( - G_ik * np.sin(theta_ik) + B_ik * np.cos(theta_ik) )
            else:
                G_ij     = Y_bus[bus_idx_i, bus_idx_j].real
                B_ij     = Y_bus[bus_idx_i, bus_idx_j].imag
                V_i      = V[bus_idx_i, 0]
                V_j      = V[bus_idx_j, 0]
                theta_ij = V[bus_idx_i, 1] - V[bus_idx_j, 1]
                J[i, j] += V_i * V_j * ( G_ij * np.sin(theta_ij) - B_ij * np.cos(theta_ij) )

    # 2) Real, differentiated by voltage magnitude (J12)
    for i, bus_idx_i in enumerate(map_idx[:(N_pq+N_pv)]):
        for j, bus_idx_j in enumerate(map_idx[(N_pq+N_pv):]):
            col = j + (N_pq+N_pv)
            
            if bus_idx_i == bus_idx_j:
                for k in range(N_bus):
                    if k == bus_idx_i:
                        J[i, col] += 2 * V[bus_idx_i, 0] * Y_bus[bus_idx_i, bus_idx_i].real
                    else:
                        G_ik     = Y_bus[bus_idx_i, k].real
                        B_ik     = Y_bus[bus_idx_i, k].imag
                        V_k      = V[k, 0]
                        theta_ik = V[bus_idx_i, 1] - V[k, 1]
                        J[i, col] += V_k * ( G_ik * np.cos(theta_ik) + B_ik * np.sin(theta_ik) )
            else:
                G_ij     = Y_bus[bus_idx_i, bus_idx_j].real
                B_ij     = Y_bus[bus_idx_i, bus_idx_j].imag
                V_i      = V[bus_idx_i, 0]
                theta_ij = V[bus_idx_i, 1] - V[bus_idx_j, 1]
                
                J[i, col] += V_i * ( G_ij * np.cos(theta_ij) + B_ij * np.sin(theta_ij) )

    # 3) Imag, differentiated by voltage angle (J21)
    for i, bus_idx_i in enumerate(map_idx[(N_pq+N_pv):]):
        row = i + (N_pq+N_pv)
        for j, bus_idx_j in enumerate(map_idx[:(N_pq+N_pv)]):
            if bus_idx_i == bus_idx_j:
                for k in range(N_bus):
                    if k == bus_idx_i: continue
                    G_ik     = Y_bus[bus_idx_i, k].real
                    B_ik     = Y_bus[bus_idx_i, k].imag
                    V_i      = V[bus_idx_i, 0]
                    V_k      = V[k, 0]
                    theta_ik = V[bus_idx_i, 1] - V[k, 1]
                    J[row, j] += V_i * V_k * ( G_ik * np.cos(theta_ik) + B_ik * np.sin(theta_ik) )
            else:
                G_ij     = Y_bus[bus_idx_i, bus_idx_j].real
                B_ij     = Y_bus[bus_idx_i, bus_idx_j].imag
                V_i      = V[bus_idx_i, 0]
                V_j      = V[bus_idx_j, 0]
                theta_ij = V[bus_idx_i, 1] - V[bus_idx_j, 1]
                J[row, j] += V_i * V_j * ( - G_ij * np.cos(theta_ij) - B_ij * np.sin(theta_ij) )

    # 4) Imag, differentiated by voltage magnitude (J22)
    for i, bus_idx_i in enumerate(map_idx[(N_pq+N_pv):]):
        row = i + (N_pq+N_pv)
        for j, bus_idx_j in enumerate(map_idx[(N_pq+N_pv):]):
            col = j + (N_pq+N_pv)
            
            if bus_idx_i == bus_idx_j:
                for k in range(N_bus):
                    if k == bus_idx_i:
                        J[row, col] -= 2 * V[bus_idx_i, 0] * Y_bus[bus_idx_i, bus_idx_i].imag
                    else:
                        G_ik     = Y_bus[bus_idx_i, k].real
                        B_ik     = Y_bus[bus_idx_i, k].imag
                        V_k      = V[k, 0]
                        theta_ik = V[bus_idx_i, 1] - V[k, 1]
                        J[row, col] += V_k * ( G_ik * np.sin(theta_ik) - B_ik * np.cos(theta_ik) )
            else:
                G_ij     = Y_bus[bus_idx_i, bus_idx_j].real
                B_ij     = Y_bus[bus_idx_i, bus_idx_j].imag
                V_i      = V[bus_idx_i, 0]
                theta_ij = V[bus_idx_i, 1] - V[bus_idx_j, 1]

                J[row, col] += V_i * ( G_ij * np.sin(theta_ij) - B_ij * np.cos(theta_ij) )
    
    return J



def Mismatch( map_tuple, Y_bus, V, Bus_Data, Gen_Data, N_bus, dim, N_pq, N_pv, S_base ):

    N_gen = Gen_Data.shape[0]
    mismatch_total = np.zeros( (N_bus, 2) )

    f_x = np.zeros( dim )

    for i in range( N_bus ):

        mismatch_total[i, 0] += Bus_Data[i, 2] / S_base  # Load P (MW -> p.u.)
        mismatch_total[i, 1] += Bus_Data[i, 3] / S_base  # Load Q (MVAR -> p.u.)

        for k in range( N_bus ):
            G_ik     = Y_bus[i, k].real
            B_ik     = Y_bus[i, k].imag
            V_i      = V[i, 0]
            V_k      = V[k, 0]
            theta_ik = V[i, 1] - V[k, 1]

            mismatch_total[i, 0]  += V_i * V_k * ( G_ik * np.cos(theta_ik) + B_ik * np.sin(theta_ik) )   # Real
            mismatch_total[i, 1]  += V_i * V_k * ( G_ik * np.sin(theta_ik) - B_ik * np.cos(theta_ik) )   # Reactive

    for i in range( N_gen ):
        bus_idx = int( Gen_Data[i, 0] - 1 )

        mismatch_total[bus_idx, 0] -= Gen_Data[i, 1] / S_base   # Gen P (MW -> p.u.)
        mismatch_total[bus_idx, 1] -= Gen_Data[i, 2] / S_base   # Gen Q (MVAR -> p.u.)

    map_arr = np.array(map_tuple, dtype=int)

    f_x[:(N_pq+N_pv)] = mismatch_total[ map_arr[:(N_pq+N_pv)] - 1 , 0]
    f_x[(N_pq+N_pv):] = mismatch_total[ map_arr[(N_pq+N_pv):] - 1 , 1]

    return f_x



def Newton_Raphson(x_prev, J_prev, f_prev):
    # Dimension Check
    if J_prev.shape[0] != J_prev.shape[1]:
        raise ValueError(f"Jacobian is not square-form: {J_prev.shape}")
    if J_prev.shape[1] != f_prev.shape[0]:
        raise ValueError(f"Dimension mismatch: J={J_prev.shape}, f={f_prev.shape}")
    if x_prev.shape[0] != f_prev.shape[0]:
        raise ValueError(f"Dimension mismatch: x={x_prev.shape}, f={f_prev.shape}")
    
    try:
        delta_x = - np.linalg.solve(J_prev, f_prev)
    except np.linalg.LinAlgError:
        raise ValueError("Jacobian Matrix is singular/ill-conditioned.")
        
    x_next  = x_prev + delta_x
 
    return x_next, delta_x