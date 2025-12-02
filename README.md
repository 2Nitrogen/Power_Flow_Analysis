# Power Flow Analysis Using Newton-Raphson Method

This repository provides how power flow analysis works via numerical iterative method, particularly Newton-Raphson (NR) method.

This program and its code repository were developed as part of the term project for *Introduction to Electric Power and Energy Systems* at Seoul National University (2nd semester, 2025) and are being released publicly.

---
### Input Data Format

Input data files consist of three different `.csv` files, including bus, transmission line (branch), and generator data.

A set of the three data is defined for each case. Here, we use the following directory and file naming convention for case 'X'.

- CASE_X_INPUT/CASE_X_BUS DATA.csv
- CASE_X_INPUT/CASE_X_BRANCH DATA.csv
- CASE_X_INPUT/CASE_X_GEN DATA.csv

This code is written assuming that all data follow `matpower`'s custom.

Referring to this [<u>manual</u>](https://matpower.org/docs/ref/matpower7.1/lib/idx_bus.html) would be helpful.

---
### Power Flow Analysis

This sections offers detailed explanation on calculations involved in power flow analysis.

1. Power Balance Equation

$P_{G_i}-P_{D_i} = \sum_{i=1}^{N}{|V_i||V_k|\left(G_s \cos{\theta_{ik} + B_s \sin{\theta_{ik}}}\right)}$
    
$Q_{G_i}-Q_{D_i} = \sum_{i=1}^{N}{|V_i||V_k|(G_s \sin{\theta_{ik} - B_s \cos{\theta_{ik}}})}$

2. Admittance Matrix

: Based on Lumped-Equivalent Circuit Model, each element of admittance matrix $\mathbf{Y}_{bus} \in \mathbb{C}^{N}$ can be computed as follows

$$
\left( \mathbf{Y}_{\text{bus}} \right)_{ii}
= \sum_{\substack{k=1 \\ k \neq i}}^{N}
  \left( \frac{1}{z_{\text{series}}^{ik}} + \frac{y_{\text{sh}}^{ik}}{2} \right)
$$

$$
\left( \mathbf{Y}_{\text{bus}} \right)_{ij}
= - \frac{1}{z_{\text{series}}^{ij}} \quad \text{for } i \neq j
$$

: In addition, we need to consider shunt power at each bus, specified as Gs and Bs at bus data of .csv file. Each value is divided by $V_{base}^2$ to convert to admittance. Sign is determined by the custom.

$ \left(\mathbf{Y}_{bus} \right)_{ii} -= \dfrac{G_s + jB_s}{V_{base}^2} $

3. Definition of Variable for NR

: Note that numbering on each variable element may not correspond to the bus number. That's because we define voltage angle $\theta$ as variables only for PQ and PV buses, and voltage magnitude $|V|$ as variables only for PQ buses for computational efficiency. This also  applies to mismatch vector $\mathbf{f(x)}$ - numbering on power may NOT refer to bus number.

$$
\mathbf{x} =
\begin{bmatrix}
  \theta_1 \\
  \vdots \\
  \theta_{N_{pq} + N_{pv}} \\
  |V|_1 \\
  \vdots \\
  |V|_{N_{pq}}
\end{bmatrix}
\in \mathbb{R}^{2 N_{pq} + N_{pv}}
$$

$$
\mathbf{f}(\mathbf{x}) =
\begin{bmatrix}
  P_1 - P_{G_1} + P_{D_1} \\
  \vdots \\
  P_{N_{pq} + N_{pv}} - P_{G_{N_{pq} + N_{pv}}} + P_{D_{N_{pq} + N_{pv}}} \\
  Q_1 - Q_{G_1} + Q_{D_1} \\
  \vdots \\
  Q_{N_{pq}} - Q_{G_{N_{pq}}} + Q_{D_{N_{pq}}}
\end{bmatrix}
\in \mathbb{R}^{2 N_{pq} + N_{pv}}
$$


4. Jacobian

: we need to calculate jacobian of mismatch vector $\mathbf{f(x)}$ differentiated by variable $\mathbf{x}$ for NR iteration. We can divide the jacobian into four different sections as below.

$\mathbf{J(x)} := \dfrac{\partial \mathbf{f(x)}}{\partial \mathbf{x}} = \begin{bmatrix} \frac{\partial P}{\partial \theta}, \frac{\partial P}{\partial |V|} \\ \frac{\partial Q}{\partial \theta}, \frac{\partial Q}{\partial |V|} \end{bmatrix} = \begin{bmatrix} \mathbf{J}_{11}, \mathbf{J}_{12} \\ \mathbf{J}_{21}, \mathbf{J}_{22} \end{bmatrix}$

4-a. Real Power differentiated by voltage angle ($\frac{\partial P}{\partial \theta}$, $\mathbf{J}_{11}$)

$\frac{\partial P_i}{\partial \theta_i} = \sum_{k=1, k \neq i}^{N}{|V_i||V_k|(-G_{ik} \sin{\theta_{ik} + B_{ik} \cos{\theta_{ik}}})} = -Q_i-B_{ii}|V_i|^2$

$\frac{\partial P_i}{\partial \theta_j} = |V_i||V_j|(G_{ij} \sin{\theta_{ij} - B_{ij} \cos{\theta_{ij}}}) \quad \small{(i \neq j)}$

4-b. Real Power differentiated by voltage magnitude ($\frac{\partial P}{\partial |V|}$, $\mathbf{J}_{12}$)

$\frac{\partial P_i}{\partial |V_i|} = \sum_{k=1, k \neq i}^{N}{|V_k|(G_{ik} \cos{\theta_{ik} + B_{ik} \sin{\theta_{ik}}})} + 2|V_i|G_{ii}$

$\frac{\partial P_i}{\partial |V_j|} = |V_i|(G_{ij} \cos{\theta_{ij} + B_{ij} \sin{\theta_{ij}}}) \quad \small{(i \neq j)}$

4-c. Imaginary Power differentiated by voltage angle ($\frac{\partial Q}{\partial \theta}$, $\mathbf{J}_{21}$)

$\frac{\partial Q_i}{\partial \theta_i} = \sum_{k=1, k \neq i}^{N}{|V_i||V_k|(G_{ik} \cos{\theta_{ik} + B_{ik} \sin{\theta_{ik}}})} = P_i-G_{ii}|V_i|^2$

$\frac{\partial Q_i}{\partial \theta_j} = |V_i||V_j|(-G_{ij} \cos{\theta_{ij} - B_{ij} \sin{\theta_{ij}}}) \quad \small{(i \neq j)}$

4-d. Imaginary Power differentiated by voltage magnitude ($\frac{\partial Q}{\partial |V|}$, $\mathbf{J}_{22}$)

$\frac{\partial Q_i}{\partial |V_i|} = \sum_{k=1, k \neq i}^{N}{|V_k|(G_{ik} \sin{\theta_{ik} - B_{ik} \cos{\theta_{ik}}})} - 2|V_i|B_{ii}$

$\frac{\partial Q_i}{\partial |V_j|} = |V_i|(G_{ij} \sin{\theta_{ij} - B_{ij} \cos{\theta_{ij}}}) \quad \small{(i \neq j)}$


5. NR iteration

: implement the following iteration rule until termination criteria is satisfied, or until the total number of iteration does not exceed predefined maximum iteration number.

$$
\mathbf{x}^{(s+1)} = \mathbf{x}^{(s)}-\mathbf{J(\mathbf{x}^{(s)})}^{-1} \mathbf{f(\mathbf{x}^{(s)})}
$$

---
### Note

- To determine starting point of variable, $\mathbf{x}^{(0)}$, Fast Decoupled Power Flow method has been implemented. However, for small systems without ill-conditioning of Jacobian matrix, it is much better to use flat starting point (voltage magnitude=1p.u., voltage angle=0).

- To prevent any dependency issues, please refer to configurations below.
    * Python = 3.12.x
    * Numpy  = 2.0.2
    * Pandas = 2.2.3