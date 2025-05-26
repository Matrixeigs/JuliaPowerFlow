# Power System Component Mathematical Models

This document provides comprehensive mathematical formulations for all power system components implemented in the JuliaPowerFlow package.

## 1. Bus Models

### 1.1 Bus Types and Equations

Buses are nodes in the power system where components connect. They are classified into three types:

#### 1.1.1 PQ Bus (Load Bus)
- **Known**: Active power $P$ and reactive power $Q$
- **Unknown**: Voltage magnitude $|V|$ and angle $\theta$
- **Equations**:
$$
P_i = |V_i| \sum_{j=1}^n |V_j| |Y_{ij}| \cos(\theta_i - \theta_j - \alpha_{ij}) = P_{i,\text{specified}}
$$
$$
Q_i = |V_i| \sum_{j=1}^n |V_j| |Y_{ij}| \sin(\theta_i - \theta_j - \alpha_{ij}) = Q_{i,\text{specified}}
$$

#### 1.1.2 PV Bus (Generator Bus)
- **Known**: Active power $P$ and voltage magnitude $|V|$
- **Unknown**: Reactive power $Q$ and voltage angle $\theta$
- **Equations**:
$$
P_i = |V_i| \sum_{j=1}^n |V_j| |Y_{ij}| \cos(\theta_i - \theta_j - \alpha_{ij}) = P_{i,\text{specified}}
$$
$$
|V_i| = |V_{i,\text{specified}}|
$$

#### 1.1.3 Slack Bus (Reference Bus)
- **Known**: Voltage magnitude $|V|$ and angle $\theta$ (typically $\theta = 0$)
- **Unknown**: Active power $P$ and reactive power $Q$
- **Equations**:
$$
|V_i| = |V_{i,\text{specified}}|, \quad \theta_i = \theta_{i,\text{specified}}
$$

### 1.2 Bus Power Balance

At each bus, Kirchhoff's current law gives:
$$
I_i = \sum_{j=1}^n Y_{ij} V_j = Y_{ii} V_i + \sum_{j \neq i} Y_{ij} V_j
$$

The complex power injection is:
$$
S_i = P_i + jQ_i = V_i I_i^* = V_i \left(\sum_{j=1}^n Y_{ij}^* V_j^*\right)
$$

### 1.3 Bus Shunt Elements

For buses with shunt admittances:
$$
Y_{\text{shunt},i} = G_{\text{shunt},i} + jB_{\text{shunt},i}
$$

The shunt power consumption is:
$$
S_{\text{shunt},i} = V_i^2 Y_{\text{shunt},i}^*
$$

## 2. Branch Models

### 2.1 Transmission Line Model

#### 2.1.1 π-Model Representation
A transmission line is represented by the π-equivalent circuit:

$$
\begin{bmatrix}
I_i \\
I_j
\end{bmatrix} = 
\begin{bmatrix}
y + \frac{y_c}{2} & -y \\
-y & y + \frac{y_c}{2}
\end{bmatrix}
\begin{bmatrix}
V_i \\
V_j
\end{bmatrix}
$$

where:
- $y = \frac{1}{r + jx}$ is the series admittance
- $y_c = jb$ is the shunt admittance
- $r$ is resistance, $x$ is reactance, $b$ is susceptance

#### 2.1.2 Branch Power Flows
Power flow from bus $i$ to bus $j$:
$$
S_{ij} = V_i I_{ij}^* = V_i \left[y^*(V_i^* - V_j^*) + \frac{y_c^*}{2}V_i^*\right]
$$

Power flow from bus $j$ to bus $i$:
$$
S_{ji} = V_j I_{ji}^* = V_j \left[y^*(V_j^* - V_i^*) + \frac{y_c^*}{2}V_j^*\right]
$$

#### 2.1.3 Line Losses
Total losses in the branch:
$$
S_{\text{loss}} = S_{ij} + S_{ji} = |I_{ij}|^2 (r + jx)
$$

### 2.2 Transformer Model

#### 2.2.1 Ideal Transformer with Tap Changer
For a transformer with complex tap ratio $t = |t|e^{j\phi}$:

$$
\begin{bmatrix}
I_i \\
I_j
\end{bmatrix} = 
\begin{bmatrix}
\frac{y}{|t|^2} & -\frac{y}{t^*} \\
-\frac{y}{t} & y
\end{bmatrix}
\begin{bmatrix}
V_i \\
V_j
\end{bmatrix}
$$

where:
- $|t|$ is the tap magnitude
- $\phi$ is the phase shift angle

#### 2.2.2 Transformer Power Flows
Power flow equations become:
$$
S_{ij} = V_i \left[\frac{y^*}{|t|^2}V_i^* - \frac{y^*}{t}V_j^*\right]
$$
$$
S_{ji} = V_j \left[y^*V_j^* - \frac{y^*}{t^*}V_i^*\right]
$$

#### 2.2.3 OLTC (On-Load Tap Changer) Control
The tap ratio is adjusted to maintain voltage:
$$
|V_j| = |V_{j,\text{target}}| \Rightarrow |t| = \frac{|V_i|}{|V_{j,\text{target}}|}
$$

## 3. Generator Models

### 3.1 Static Generator Model

#### 3.1.1 Basic Power Equations
For a cylindrical rotor synchronous generator:
$$
P = \frac{E_q U_t}{X_d} \sin \delta
$$
$$
Q = \frac{E_q U_t \cos \delta - U_t^2}{X_d}
$$

where:
- $E_q$ is internal voltage magnitude
- $U_t$ is terminal voltage magnitude  
- $X_d$ is direct-axis reactance
- $\delta$ is power angle

#### 3.1.2 Capability Constraints
The generator operating region is bounded by:

**Apparent Power Limit:**
$$
P^2 + Q^2 \leq S_{\max}^2
$$

**Active Power Limits:**
$$
P_{\min} \leq P \leq P_{\max}
$$

**Reactive Power Limits:**
$$
Q_{\min} \leq Q \leq Q_{\max}
$$

**Power Angle Limits:**
$$
|\delta| \leq \delta_{\max}
$$

**Field Current Limits:**
$$
E_{q,\min} \leq E_q \leq E_{q,\max}
$$

#### 3.1.3 Dynamic Q_min Function
The minimum reactive power varies with active power:
$$
Q_{\min}(P) = \max\left(Q_{\min,\text{constant}}, \frac{E_{q,\min}U_t\cos(\delta_{\max}) - U_t^2}{X_d}\right)
$$

### 3.2 Generator Control Models

#### 3.2.1 Automatic Voltage Regulator (AVR)
$$
E_q = E_{q0} + K_A(V_{\text{ref}} - V_t)
$$

#### 3.2.2 Governor Control
$$
P = P_0 + \frac{1}{R}(f_{\text{ref}} - f)
$$

where $R$ is the speed regulation (droop).

## 4. Load Models

### 4.1 Static Load Models

#### 4.1.1 Constant Power Load
$$
P = P_0, \quad Q = Q_0
$$

#### 4.1.2 Constant Current Load
$$
P = P_0 \frac{V}{V_0}, \quad Q = Q_0 \frac{V}{V_0}
$$

#### 4.1.3 Constant Impedance Load
$$
P = P_0 \left(\frac{V}{V_0}\right)^2, \quad Q = Q_0 \left(\frac{V}{V_0}\right)^2
$$

#### 4.1.4 ZIP Load Model
Combination of constant impedance (Z), current (I), and power (P):
$$
P = P_0 \left[a_1 \left(\frac{V}{V_0}\right)^2 + a_2 \frac{V}{V_0} + a_3\right]
$$
$$
Q = Q_0 \left[b_1 \left(\frac{V}{V_0}\right)^2 + b_2 \frac{V}{V_0} + b_3\right]
$$

where $a_1 + a_2 + a_3 = 1$ and $b_1 + b_2 + b_3 = 1$.

#### 4.1.5 Exponential Load Model
$$
P = P_0 \left(\frac{V}{V_0}\right)^{\alpha_p}
$$
$$
Q = Q_0 \left(\frac{V}{V_0}\right)^{\alpha_q}
$$

### 4.2 Motor Load Models

#### 4.2.1 Induction Motor Equivalent Circuit
$$
P_{\text{motor}} = \frac{3V^2 r_2/s}{(r_1 + r_2/s)^2 + (x_1 + x_2)^2}
$$

where $s$ is the slip and $r_1, x_1, r_2, x_2$ are circuit parameters.

## 5. System-Level Models

### 5.1 Y-bus Formation

#### 5.1.1 Nodal Admittance Matrix
The bus admittance matrix relates bus currents to voltages:
$$
\mathbf{I} = \mathbf{Y}_{\text{bus}} \mathbf{V}
$$

#### 5.1.2 Y-bus Construction Algorithm
1. **Initialize**: $Y_{ij} = 0$ for all $i,j$
2. **Add branch admittances**:
   - $Y_{ii} = Y_{ii} + y_{ij} + \frac{y_{c,ij}}{2}$
   - $Y_{jj} = Y_{jj} + y_{ij} + \frac{y_{c,ij}}{2}$
   - $Y_{ij} = Y_{ij} - y_{ij}$
   - $Y_{ji} = Y_{ji} - y_{ij}$
3. **Add shunt elements**: $Y_{ii} = Y_{ii} + y_{\text{shunt},i}$

### 5.2 Power Flow Equations

#### 5.2.1 Complex Power Balance
At each bus $i$:
$$
S_i = P_i + jQ_i = V_i \sum_{j=1}^n Y_{ij}^* V_j^*
$$

#### 5.2.2 Rectangular Form
In rectangular coordinates $V_i = e_i + jf_i$:
$$
P_i = e_i \sum_{j=1}^n (G_{ij}e_j - B_{ij}f_j) + f_i \sum_{j=1}^n (G_{ij}f_j + B_{ij}e_j)
$$
$$
Q_i = f_i \sum_{j=1}^n (G_{ij}e_j - B_{ij}f_j) - e_i \sum_{j=1}^n (G_{ij}f_j + B_{ij}e_j)
$$

#### 5.2.3 Polar Form
In polar coordinates $V_i = |V_i|e^{j\theta_i}$:
$$
P_i = |V_i| \sum_{j=1}^n |V_j| |Y_{ij}| \cos(\theta_i - \theta_j - \alpha_{ij})
$$
$$
Q_i = |V_i| \sum_{j=1}^n |V_j| |Y_{ij}| \sin(\theta_i - \theta_j - \alpha_{ij})
$$

### 5.3 Jacobian Matrix

#### 5.3.1 Newton-Raphson Jacobian
The Jacobian matrix for power flow has the structure:
$$
\mathbf{J} = \begin{bmatrix}
\frac{\partial \mathbf{P}}{\partial \boldsymbol{\theta}} & \frac{\partial \mathbf{P}}{\partial |\mathbf{V}|} \\
\frac{\partial \mathbf{Q}}{\partial \boldsymbol{\theta}} & \frac{\partial \mathbf{Q}}{\partial |\mathbf{V}|}
\end{bmatrix}
$$

#### 5.3.2 Jacobian Elements
$$
\frac{\partial P_i}{\partial \theta_j} = \begin{cases}
|V_i| |V_j| |Y_{ij}| \sin(\theta_i - \theta_j - \alpha_{ij}) & i \neq j \\
-Q_i - |V_i|^2 B_{ii} & i = j
\end{cases}
$$

$$
\frac{\partial P_i}{\partial |V_j|} = \begin{cases}
|V_i| |Y_{ij}| \cos(\theta_i - \theta_j - \alpha_{ij}) & i \neq j \\
\frac{P_i}{|V_i|} + |V_i| G_{ii} & i = j
\end{cases}
$$

$$
\frac{\partial Q_i}{\partial \theta_j} = \begin{cases}
-|V_i| |V_j| |Y_{ij}| \cos(\theta_i - \theta_j - \alpha_{ij}) & i \neq j \\
P_i - |V_i|^2 G_{ii} & i = j
\end{cases}
$$

$$
\frac{\partial Q_i}{\partial |V_j|} = \begin{cases}
|V_i| |Y_{ij}| \sin(\theta_i - \theta_j - \alpha_{ij}) & i \neq j \\
\frac{Q_i}{|V_i|} - |V_i| B_{ii} & i = j
\end{cases}
$$

## 6. Control System Models

### 6.1 Frequency Control

#### 6.1.1 Primary Control (Governor)
$$
\Delta P_m = -\frac{1}{R} \Delta f
$$

#### 6.1.2 Secondary Control (AGC)
$$
\Delta P_{\text{ref}} = -K_I \int \text{ACE} \, dt
$$

where ACE is the Area Control Error.

### 6.2 Voltage Control

#### 6.2.1 Automatic Voltage Regulator
$$
\frac{dE_{fd}}{dt} = \frac{1}{T_A}(K_A(V_{\text{ref}} - V_t) - E_{fd})
$$

#### 6.2.2 Tap Changer Control
$$
\frac{dn}{dt} = \frac{1}{T_{\text{delay}}} \text{sign}(V_{\text{target}} - V_{\text{controlled}})
$$

## 7. Implementation Considerations

### 7.1 Numerical Stability
- Use appropriate per-unit system
- Handle singular cases in matrix operations
- Implement proper convergence criteria

### 7.2 Computational Efficiency
- Exploit sparsity in Y-bus matrix
- Use factorization techniques for repeated solutions
- Implement warm-start capabilities

### 7.3 Model Validation
- Compare with analytical solutions for simple cases
- Validate against commercial software
- Cross-check with field measurements

This comprehensive mathematical framework provides the foundation for all component models implemented in the JuliaPowerFlow package.
