import numpy as np   # 2.0.2
from tools import *

##### 1. Load Data #####
# [NOTE] Change this factor to check other cases (ex. "A", "B", or "C")
CASE = "A"
Bus_Data, Branch_Data, Gen_Data = Read_Data(CASE)

##### 2. Define Variables & Jacobian #####
# Dimensions
N_bus = Bus_Data.shape[0]
N_pq = np.sum(Bus_Data[:, 1] == 1)
N_pv = np.sum(Bus_Data[:, 1] == 2)
N_slack = np.sum(Bus_Data[:, 1] == 3)

if N_bus != (N_pq + N_pv + N_slack):
    raise ValueError("Given Bus Data doesn't fit the standard data format.")


V = np.zeros((N_bus, 2))  # following bus number order
gen_v_map = dict(zip(Gen_Data[:, 0], Gen_Data[:, 3]))

# Fill up determined Voltage using data sheet
for i in range(N_bus):
    bus_num = int(Bus_Data[i, 0])
    bus_type = Bus_Data[i, 1]  # 1=PQ, 2=PV, 3=Slack

    if bus_type == 3:    # Slack Bus
        if bus_num in gen_v_map:
            V[i, 0] = gen_v_map[bus_num]  # Slack voltage magnitude
        else:
            V[i, 0] = 1.0
        V[i, 1] = 0.0

    elif bus_type == 2:  # PV Bus
        if bus_num in gen_v_map:
            V[i, 0] = gen_v_map[bus_num]
        else:
            raise ValueError(f"Bus {bus_num} is PV but no Gen Data found.")

    elif bus_type == 1:  # PQ Bus
        continue

    else:
        raise ValueError(f"Invalid bus type at row {i}: {bus_type}")


dim = 2 * N_pq + N_pv
x = np.zeros(dim)  # consider non-slack buses only

temp = [None] * dim

pq_mask = Bus_Data[:, 1] == 1
pv_mask = Bus_Data[:, 1] == 2

temp[:N_pq] = Bus_Data[pq_mask, 0]
temp[N_pq:N_pq + N_pv] = Bus_Data[pv_mask, 0]
temp[N_pq + N_pv:] = Bus_Data[pq_mask, 0]
t = tuple(temp)   # Mapping tuple for variable vector 'x' ... ( PQ - PV )


S_base = 100  # Base power for per-unit analysis [MVA]

Y_bus = admittance_matrix(
    N_bus,
    Bus_Data,
    Branch_Data,
    S_base,
)  # Admittance Matrix

J = np.zeros((dim, dim))  # Jacobian
f = np.zeros(dim)         # mismatch vector


##### 3. Power Flow Analysis #####

max_iter = 1000  # max number of iteration
eps = 1e-4       # termination criteria

# Initial guess using FDPF (iterations = 0: Flat Start)
# [NOTE] Here, because our system is very small, flat start is much better
# (may consider using FDPF for systems with ill-conditioned Jacobian!)
x0 = np.zeros_like(x)
x0 = Approx_starting_point(
    Y_bus=Y_bus,
    Bus_Data=Bus_Data,
    Gen_Data=Gen_Data,
    S_base=S_base,
    N_bus=N_bus,
    N_pq=N_pq,
    N_pv=N_pv,
    map_tuple=t,
    dim=dim,
    iterations=0,
)

# NR method
for i in range(max_iter):

    for p in range(N_pq + N_pv):
        V[t[p] - 1, 1] = x0[p]

    for q in range(N_pq):
        V[t[q + (N_pq + N_pv)] - 1, 0] = x0[q + (N_pq + N_pv)]

    J = Jacobian(
        map_tuple=t,
        Y_bus=Y_bus,
        V=V,
        N_bus=N_bus,
        dim=dim,
        N_pq=N_pq,
        N_pv=N_pv,
    )

    f = Mismatch(
        map_tuple=t,
        Y_bus=Y_bus,
        V=V,
        Bus_Data=Bus_Data,
        Gen_Data=Gen_Data,
        N_bus=N_bus,
        dim=dim,
        N_pq=N_pq,
        N_pv=N_pv,
        S_base=S_base,
    )

    x, delta_x = Newton_Raphson(x0, J_prev=J, f_prev=f)
    x0 = x

    if np.max(np.abs(delta_x)) < eps:
        print(f"Termination Criteria satisfied ... Total Number of Iteration = {i + 1}")
        break


##### 5. Print Results #####
for i in range(N_bus):
    print(
        f"Bus No. {i + 1}: "
        f"Magnitude = {V[i, 0]:.4f} [p.u.], "
        f"Angle = {V[i, 1] * 180 / np.pi:.4f} [deg]"
    )
