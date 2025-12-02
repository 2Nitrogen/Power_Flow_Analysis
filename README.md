# Power Flow Analysis Using Newton-Raphson Method

This repository demonstrates how power flow analysis works via numerical iterative methods, particularly the Newton-Raphson (NR) method.

This program and its code repository were developed as part of the term project for *Introduction to Electric Power and Energy Systems* at Seoul National University (2nd semester, 2025) and are being released publicly.

---

### Input Data Format

Input data files consist of three different `.csv` files: bus data, transmission line (branch) data, and generator data.

A set of these three files is defined for each case. For a case `X`, we use the following directory and file naming convention:

- `CASE_X_INPUT/CASE_X_BUS.csv`
- `CASE_X_INPUT/CASE_X_BRANCH.csv`
- `CASE_X_INPUT/CASE_X_GEN.csv`

This code is written assuming that all data follow **MATPOWER**’s convention.

Referring to this [manual](https://matpower.org/docs/ref/matpower7.1/lib/idx_bus.html) would be helpful.

---

### Power Flow Analysis

This section offers a detailed explanation of the calculations involved in power flow analysis.

---

#### 1. Power Balance Equations

For each bus \(i\), the active and reactive power balance equations are

$$
P_{G_i} - P_{D_i}
= \sum_{k=1}^{N} |V_i||V_k|
\left( G_{ik} \cos \theta_{ik} + B_{ik} \sin \theta_{ik} \right),
$$

$$
Q_{G_i} - Q_{D_i}
= \sum_{k=1}^{N} |V_i||V_k|
\left( G_{ik} \sin \theta_{ik} - B_{ik} \cos \theta_{ik} \right),
$$

where \(G_{ik} + jB_{ik}\) is the \((i,k)\)-th element of the bus admittance matrix \(\mathbf{Y}_{\text{bus}}\), and \(\theta_{ik} = \theta_i - \theta_k\).

---

#### 2. Admittance Matrix

Based on the lumped-equivalent \(\pi\)-model, each element of the bus admittance matrix
\(\mathbf{Y}_{\text{bus}} \in \mathbb{C}^{N \times N}\) can be computed as follows.

For a transmission line between bus \(i\) and bus \(k\) with series impedance \(z_{\text{series}}^{ik}\) and total shunt admittance \(y_{\text{sh}}^{ik}\),

$$
(\mathbf{Y}_{\text{bus}})_{ii}
= \sum_{\substack{k=1 \\ k \neq i}}^{N}
\left(
\frac{1}{z_{\text{series}}^{ik}} + \frac{y_{\text{sh}}^{ik}}{2}
\right),
$$

$$
(\mathbf{Y}_{\text{bus}})_{ij}
= -\,\frac{1}{z_{\text{series}}^{ij}},
\qquad i \neq j.
$$

In addition, we need to consider shunt power at each bus, specified as \(G_s\) and \(B_s\) in the bus data `.csv` file. Each value is divided by \(V_{\text{base}}^2\) to convert to admittance. The sign is determined by the adopted convention.

$$
(\mathbf{Y}_{\text{bus}})_{ii}
\;\mathrel{-}= \frac{G_s + j B_s}{V_{\text{base}}^2}.
$$

---

#### 3. Definition of Variables for Newton-Raphson

Note that the ordering of the elements in the state vector does **not** necessarily correspond to the bus numbers. This is because we define

- voltage angles \(\theta\) as variables only for PQ and PV buses, and  
- voltage magnitudes \(|V|\) as variables only for PQ buses,

for computational efficiency.

The same reordering applies to the mismatch vector \(\mathbf{f}(\mathbf{x})\); indices in \(\mathbf{f}(\mathbf{x})\) do **not** directly represent bus numbers.

The state vector and mismatch vector are defined as

$$
\mathbf{x}
=
\begin{bmatrix}
\theta_1 \\
\vdots \\
\theta_{N_{pq} + N_{pv}} \\
|V|_1 \\
\vdots \\
|V|_{N_{pq}}
\end{bmatrix}
\in \mathbb{R}^{\,2N_{pq} + N_{pv}},
$$

$$
\mathbf{f}(\mathbf{x})
=
\begin{bmatrix}
P_1 - P_{G_1} + P_{D_1} \\
\vdots \\
P_{N_{pq} + N_{pv}} - P_{G_{N_{pq} + N_{pv}}} + P_{D_{N_{pq} + N_{pv}}} \\
Q_1 - Q_{G_1} + Q_{D_1} \\
\vdots \\
Q_{N_{pq}} - Q_{G_{N_{pq}}} + Q_{D_{N_{pq}}}
\end{bmatrix}
\in \mathbb{R}^{\,2N_{pq} + N_{pv}}.
$$

---

#### 4. Jacobian Matrix

We need to calculate the Jacobian of the mismatch vector \(\mathbf{f}(\mathbf{x})\) with respect to the state vector \(\mathbf{x}\) for the NR iteration. The Jacobian can be partitioned into four submatrices as

$$
\mathbf{J}(\mathbf{x})
:= \frac{\partial \mathbf{f}(\mathbf{x})}{\partial \mathbf{x}}
=
\begin{bmatrix}
\frac{\partial P}{\partial \theta} & \frac{\partial P}{\partial |V|} \\
\frac{\partial Q}{\partial \theta} & \frac{\partial Q}{\partial |V|}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{J}_{11} & \mathbf{J}_{12} \\
\mathbf{J}_{21} & \mathbf{J}_{22}
\end{bmatrix}.
$$

Below, \(G_{ik} + jB_{ik}\) denotes elements of \(\mathbf{Y}_{\text{bus}}\), and \(\theta_{ik} = \theta_i - \theta_k\).

---

##### 4-a. Real Power w.r.t. Voltage Angle \(\left(\dfrac{\partial P}{\partial \theta},\, \mathbf{J}_{11}\right)\)

For the diagonal terms:

$$
\frac{\partial P_i}{\partial \theta_i}
=
\sum_{\substack{k=1 \\ k \neq i}}^{N}
|V_i||V_k|
\left(
- G_{ik} \sin \theta_{ik}
+ B_{ik} \cos \theta_{ik}
\right)
=
- Q_i - B_{ii} |V_i|^2.
$$

For the off-diagonal terms \((i \neq j)\):

$$
\frac{\partial P_i}{\partial \theta_j}
=
|V_i||V_j|
\left(
G_{ij} \sin \theta_{ij}
- B_{ij} \cos \theta_{ij}
\right),
\qquad i \neq j.
$$

---

##### 4-b. Real Power w.r.t. Voltage Magnitude \(\left(\dfrac{\partial P}{\partial |V|},\, \mathbf{J}_{12}\right)\)

For the diagonal terms:

$$
\frac{\partial P_i}{\partial |V_i|}
=
\sum_{\substack{k=1 \\ k \neq i}}^{N}
|V_k|
\left(
G_{ik} \cos \theta_{ik}
+ B_{ik} \sin \theta_{ik}
\right)
+ 2 |V_i| G_{ii}.
$$

For the off-diagonal terms \((i \neq j)\):

$$
\frac{\partial P_i}{\partial |V_j|}
=
|V_i|
\left(
G_{ij} \cos \theta_{ij}
+ B_{ij} \sin \theta_{ij}
\right),
\qquad i \neq j.
$$

---

##### 4-c. Reactive Power w.r.t. Voltage Angle \(\left(\dfrac{\partial Q}{\partial \theta},\, \mathbf{J}_{21}\right)\)

For the diagonal terms:

$$
\frac{\partial Q_i}{\partial \theta_i}
=
\sum_{\substack{k=1 \\ k \neq i}}^{N}
|V_i||V_k|
\left(
G_{ik} \cos \theta_{ik}
+ B_{ik} \sin \theta_{ik}
\right)
=
P_i - G_{ii} |V_i|^2.
$$

For the off-diagonal terms \((i \neq j)\):

$$
\frac{\partial Q_i}{\partial \theta_j}
=
|V_i||V_j|
\left(
- G_{ij} \cos \theta_{ij}
- B_{ij} \sin \theta_{ij}
\right),
\qquad i \neq j.
$$

---

##### 4-d. Reactive Power w.r.t. Voltage Magnitude \(\left(\dfrac{\partial Q}{\partial |V|},\, \mathbf{J}_{22}\right)\)

For the diagonal terms:

$$
\frac{\partial Q_i}{\partial |V_i|}
=
\sum_{\substack{k=1 \\ k \neq i}}^{N}
|V_k|
\left(
G_{ik} \sin \theta_{ik}
- B_{ik} \cos \theta_{ik}
\right)
- 2 |V_i| B_{ii}.
$$

For the off-diagonal terms \((i \neq j)\):

$$
\frac{\partial Q_i}{\partial |V_j|}
=
|V_i|
\left(
G_{ij} \sin \theta_{ij}
- B_{ij} \cos \theta_{ij}
\right),
\qquad i \neq j.
$$

---

#### 5. Newton-Raphson Iteration

We implement the following iteration rule until the termination criteria are satisfied, or until the total number of iterations exceeds a predefined maximum:

$$
\mathbf{x}^{(s+1)}
=
\mathbf{x}^{(s)}
-
\mathbf{J}\big(\mathbf{x}^{(s)}\big)^{-1}
\mathbf{f}\big(\mathbf{x}^{(s)}\big),
$$

where \(s\) is the iteration index.

---

### Note

- To determine the starting point of the state vector \(\mathbf{x}^{(0)}\), a Fast Decoupled Power Flow method has been implemented. However, for small systems without severe ill-conditioning of the Jacobian matrix, it is often preferable to use a flat start
  (voltage magnitude \(= 1\) p.u., voltage angle \(= 0\)).

- To prevent any dependency issues, please refer to the configurations below:

  - Python = 3.12.x  
  - NumPy  = 2.0.2  
  - Pandas = 2.2.3
