# Mathematical Models for Electromagnetic Transient Simulation

This document outlines the mathematical models used for the components and network in the electromagnetic transient (EMT) simulation of an AC distribution system with inverter-based resources (IBRs).

## 1. Nomenclature

| Symbol | Description |
| :--- | :--- |
| $v, i$ | Instantaneous voltage and current |
| $V, I$ | Phasor or RMS voltage and current |
| $\omega$ | Angular frequency (rad/s) |
| $\delta$ | Rotor angle or power angle (rad) |
| $\psi$ | Flux linkage (Wb) |
| $L$ | Inductance (H) |
| $R$ | Resistance ($\Omega$) |
| $P, Q$ | Active and reactive power (W, VAR) |
| $H$ | Inertia constant (s) |
| $T_m, T_e$ | Mechanical and electrical torque (Nm) |
| $d, q, 0$ | Subscripts for direct, quadrature, and zero sequence components |

## 2. Reference Frame and Transformations

The simulation uses the arbitrary rotating reference frame (dq0 frame) for machine and inverter models. The network is solved in the stationary abc frame.

### Park's Transformation (abc to dq0)

The transformation from the stationary three-phase `abc` frame to a rotating `dq0` frame is given by:

$$
\begin{bmatrix} f_d \\ f_q \\ f_0 \end{bmatrix} = \frac{2}{3} \begin{bmatrix} \cos(\theta) & \cos(\theta - \frac{2\pi}{3}) & \cos(\theta + \frac{2\pi}{3}) \\ -\sin(\theta) & -\sin(\theta - \frac{2\pi}{3}) & -\sin(\theta + \frac{2\pi}{3}) \\ \frac{1}{2} & \frac{1}{2} & \frac{1}{2} \end{bmatrix} \begin{bmatrix} f_a \\ f_b \\ f_c \end{bmatrix}
$$

where $\theta = \int \omega(t) dt$ is the angle of the rotating reference frame.

### Inverse Park's Transformation (dq0 to abc)

The inverse transformation is:

$$
\begin{bmatrix} f_a \\ f_b \\ f_c \end{bmatrix} = \begin{bmatrix} \cos(\theta) & -\sin(\theta) & 1 \\ \cos(\theta - \frac{2\pi}{3}) & -\sin(\theta - \frac{2\pi}{3}) & 1 \\ \cos(\theta + \frac{2\pi}{3}) & -\sin(\theta + \frac{2\pi}{3}) & 1 \end{bmatrix} \begin{bmatrix} f_d \\ f_q \\ f_0 \end{bmatrix}
$$

## 3. Component Models

### 3.1. Synchronous Generator

A standard 6th-order model is used.

**Mechanical Equations (Swing Equation):**

$$
\frac{d\delta}{dt} = \omega_r - \omega_s
$$

$$
\frac{d\omega_r}{dt} = \frac{1}{2H} (T_m - T_e - D(\omega_r - \omega_s))
$$

where $T_e = \psi_d i_q - \psi_q i_d$.

**Electrical Equations (dq frame):**

$$
\frac{d\psi_d}{dt} = \omega_s (v_d + R_a i_d + \psi_q)
$$

$$
\frac{d\psi_q}{dt} = \omega_s (v_q + R_a i_q - \psi_d)
$$

$$
\frac{d\psi_{fd}}{dt} = \omega_s (v_{fd} - R_{fd} i_{fd})
$$

$$
\frac{d\psi_{kd}}{dt} = -\omega_s R_{kd} i_{kd}
$$

The flux linkages are related to currents by the machine inductances.

### 3.2. PV System (Inverter)

The PV system is modeled as a current-controlled voltage source inverter.

**DC Side (DC Link Capacitor):**

$$
C_{dc} \frac{dV_{dc}}{dt} = I_{pv} - I_{dc,inv}
$$

where $I_{pv}$ is the current from the PV array (often modeled as a function of irradiance and temperature) and $I_{dc,inv}$ is the DC current drawn by the inverter.

**AC Side (Inverter Control):**

The inverter injects currents ($i_d, i_q$) into the grid based on a control strategy (e.g., MPPT, voltage support). The dynamics are governed by the current controller (typically a PI controller):

$$
\frac{d\gamma_d}{dt} = i_{d,ref} - i_d
$$

$$
v_{d,inv} = K_p (i_{d,ref} - i_d) + K_i \gamma_d - \omega L_f i_q + v_d
$$

(A similar equation exists for the q-axis). $v_{d,inv}$ is the voltage synthesized by the inverter.

### 3.3. Battery Energy Storage System (BESS)

The BESS model is similar to the PV inverter, with the DC side replaced by a battery model.

**Battery Model:**

A simple internal voltage source with series resistance:

$$
V_{bat} = V_{oc}(SOC) - R_{bat} I_{bat}
$$

**State of Charge (SOC):**

$$
\frac{d(SOC)}{dt} = -\frac{I_{bat}}{Q_{max}}
$$

where $Q_{max}$ is the total capacity of the battery in Ampere-seconds.

### 3.4. Induction Motor

A 3rd-order model is common for transient studies.

**Mechanical Equation:**

$$
\frac{d\omega_r}{dt} = \frac{1}{2H} (T_e - T_{load})
$$

where $T_e = \psi_{dr} i_{qr} - \psi_{qr} i_{dr}$.

**Electrical Equations (Rotor Flux, dq frame):**

$$
\frac{d\psi_{dr}}{dt} = -\frac{R_r}{L_r} \psi_{dr} - (\omega_s - \omega_r) \psi_{qr} + \frac{R_r L_m}{L_r} i_{ds}
$$

$$
\frac{d\psi_{qr}}{dt} = -\frac{R_r}{L_r} \psi_{qr} + (\omega_s - \omega_r) \psi_{dr} + \frac{R_r L_m}{L_r} i_{qs}
$$

Stator currents ($i_{ds}, i_{qs}$) are algebraic variables determined by the network solution.

### 3.5. EV Charger

For EMT studies, an EV charger can be simplified as a constant power load or a controlled current source, similar to the inverter models but typically operating in grid-following mode.

$$
P_{ev} + jQ_{ev} = V_{bus} I_{ev}^*
$$

## 4. Network Equations

The network connects all components. It is modeled using nodal analysis.

### 4.1. Transmission Line

A nominal PI model is used for each phase:

$$
I_{send} = Y_{sh} V_{send} + \frac{1}{Z_{ser}}(V_{send} - V_{rec})
$$

$$
I_{rec} = Y_{sh} V_{rec} + \frac{1}{Z_{ser}}(V_{rec} - V_{send})
$$

where $Z_{ser}$ is the series impedance and $Y_{sh}$ is the shunt admittance.

### 4.2. Transformer

A transformer is modeled by its series leakage impedance ($Z_T$) and turns ratio ($n$).

$$
V_p = n V_s
$$

$$
I_p = \frac{1}{n} I_s
$$

The voltage drop is $V_{drop} = Z_T I_s$.

### 4.3. Nodal Analysis

The entire network is described by the nodal admittance matrix equation:

$$
\mathbf{Y_{bus} V_{bus} = I_{bus}}
$$

-   $\mathbf{Y_{bus}}$ is the bus admittance matrix, constructed from the parameters of all lines and transformers.
-   $\mathbf{V_{bus}}$ is the vector of unknown bus voltages.
-   $\mathbf{I_{bus}}$ is the vector of current injections from all generators, IBRs, and loads.

This matrix equation is solved at each time step of the simulation to find the bus voltages, which are then used to update the internal states of all connected components. For a fault, the $\mathbf{Y_{bus}}$ matrix is modified by adding a large admittance at the faulted bus.
