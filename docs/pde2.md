---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Numerical Partial Differential Equation II: Stability Analysis

+++

In the previous section, we explored the **Forward Time Centered Space (FTCS)** scheme, a fundamental explicit finite difference method for solving the linear advection equation:
\begin{align}
\frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0
\end{align}

The FTCS scheme approximates the time derivative using a forward difference and the spatial derivatives using centered differences.
\begin{align}
u_i^{n+1} = u_i^n - \frac{c \Delta t}{2 \Delta x} \left( u_{i+1}^n - u_{i-1}^n \right)
\end{align}

While straightforward to implement, the FTCS method is unconditionally unstable.
We also provided some Python code to demonstrate the application of the FTCS scheme to a sinusoidal initial condition:

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

# Parameters
c = 1.0          # Advection speed
L = 1.0          # Domain length
T = 1.25         # Total time
nx = 100         # Number of spatial points
dx = L / nx      # Spatial step size
dt = 0.001       # Initial time step size
nt = int(T / dt) # Number of time steps

# Stability parameter
sigma = c * dt / dx
print(sigma)

# Spatial grid
x = np.linspace(0, L, nx)
u_initial = np.sin(2 * np.pi * x)  # Initial condition: sinusoidal wave

# Initialize solution array
u = u_initial.copy()

# Time-stepping loop
for n in range(nt):
    # Apply periodic boundary conditions using np.roll
    u_new = u - (c * dt / (2 * dx)) * (np.roll(u, -1) - np.roll(u, 1))
    u = u_new

# Analytical solution
u_exact = np.sin(2 * np.pi * (x - c * T))

# Plotting the results
plt.figure(figsize=(12, 6))
plt.plot(x, u_initial,      label='Initial Condition', linestyle='--')
plt.plot(x, u_exact,        label='Exact Solution', linewidth=2)
plt.plot(x, u,              label='FTCS Scheme', linestyle=':', linewidth=2)
plt.xlabel('x')
plt.ylabel('u')
plt.legend()
plt.grid(True)
```

The simulation results illustrate how the FTCS scheme evolves the initial sinusoidal wave over time.
However, due to its unconditional stability, numerical instability always occur.

To overcome the limitations of the FTCS scheme, we introduce two more robust finite difference methods:
the Upwind Scheme and the Lax-Wendroff Scheme.
These methods enhance stability and accuracy, making them more suitable for solving advection-dominated problems.

+++

### 2.1 Upwind Scheme

The **Upwind Scheme** is a finite difference method specifically designed to handle advection-dominated problems more effectively than symmetric schemes like FTCS.
By incorporating the direction of wave propagation into the discretization of spatial derivatives, the upwind method enhances numerical stability and reduces non-physical oscillations.

In advection processes, information propagates in a specific direction determined by the flow velocity $c$.
The upwind scheme leverages this directional information to bias the spatial derivative approximation, ensuring that the numerical flux aligns with the physical transport direction.
This directional bias significantly improves the stability of the numerical solution, especially when dealing with sharp gradients or discontinuities.

The upwind scheme discretizes the spatial derivative based on the sign of the advection speed $c$:
*  **For $c > 0$** (flow to the right):
   \begin{align}
   \frac{\partial u}{\partial x} \approx \frac{u_i^n - u_{i-1}^n}{\Delta x}
   \end{align}
*  **For $c < 0$** (flow to the left):
   \begin{align}
   \frac{\partial u}{\partial x} \approx \frac{u_{i+1}^n - u_i^n}{\Delta x}
   \end{align}

Assuming $c > 0$ for this implementation, the **Upwind Scheme** update rule becomes:
\begin{align}
u_i^{n+1} = u_i^n - \frac{c \Delta t}{\Delta x} \left( u_i^n - u_{i-1}^n \right)
\end{align}
where:
* $u_i^n$ is the numerical approximation of $u$ at spatial index $i$ and time level $n$,
* $\Delta t$ and $\Delta x$ are the time and spatial step sizes, respectively.

The following Python code implements the upwind scheme to solve the linear advection equation for a sinusoidal initial condition.

```{code-cell} ipython3
# Parameters
c = 1.0          # Advection speed
L = 1.0          # Domain length
T = 1.25         # Total time
nx = 100         # Number of spatial points
dx = L / nx      # Spatial step size
dt = 0.001       # Initial time step size
nt = int(T / dt) # Number of time steps

# Stability parameter
sigma = c * dt / dx
print(f"Courant number (sigma): {sigma}")

# Check CFL condition
if sigma > 1:
    raise ValueError(f"CFL condition violated: sigma = {sigma} > 1. Please reduce dt or increase dx.")

# Spatial grid
x = np.linspace(0, L, nx, endpoint=False)
u_initial = np.sin(2 * np.pi * x)  # Initial condition: sinusoidal wave

# Initialize solution array
u = u_initial.copy()

# Time-stepping loop using Upwind scheme
for n in range(nt):
    # Apply periodic boundary conditions using np.roll
    u_new = u - sigma * (u - np.roll(u, 1))
    u = u_new

# Analytical solution
u_exact = np.sin(2 * np.pi * (x - c * T))

# Plotting the results
plt.figure(figsize=(12, 6))
plt.plot(x, u_initial,      label='Initial Condition', linestyle='--')
plt.plot(x, u_exact,        label='Exact Solution', linewidth=2)
plt.plot(x, u,              label='Upwind Scheme', linestyle=':', linewidth=2)
plt.xlabel('x')
plt.ylabel('u')
plt.title('Linear Advection Equation: Upwind Scheme vs Exact Solution')
plt.legend()
plt.grid(True)
```

## Non-Dimensionalization and Key Dimensionless Numbers

Non-dimensionalization is a fundamental technique in the analysis of partial differential equations (PDEs) and fluid dynamics.
By converting variables and equations into dimensionless forms, researchers and engineers can simplify complex systems, identify dominant physical effects, and generalize solutions across different scales.
This section explores the purpose of non-dimensionalization, introduces key dimensionless numbers—Reynolds Number (Re), Mach Number (Ma), and Prandtl Number (Pr)—and discusses their applications and significance in various physical systems.

+++

## Non-Dimensionalization and Key Dimensionless Numbers (15 minutes)

* Purpose of Non-Dimensionalization
  * Simplifying equations for analysis.
  * Identifying dominant physical effects in specific regimes.

* Key Dimensionless Numbers

  * Reynolds Number (Re): Ratio of inertial to viscous forces.
  * Mach Number (Ma): Compressibility effects.
  * Prandtl Number (Pr): Ratio of momentum diffusivity to thermal diffusivity.

* Applications
* How these numbers influence the behavior of physical systems.

+++

### Purpose of Non-Dimensionalization

Non-dimensionalization serves several critical purposes in the analysis and solution of PDEs:

1.  **Simplifying Equations for Analysis**

    By scaling variables and parameters, non-dimensionalization reduces the number of independent variables and parameters in the governing equations.
    This simplification often transforms complex, dimensional equations into a more manageable, dimensionless form.
    For instance, consider the non-dimensionalization of the Navier-Stokes equations:
    \begin{align}
    \frac{\partial \bar{\mathbf{u}}}{\partial \bar t} + (\bar{\mathbf{u}} \cdot \bar\nabla) \bar{\mathbf{u}} = -\bar\nabla \bar p + \frac{1}{\text{Re}} \bar\nabla^2 \bar{\mathbf{u}} + \bar{\mathbf{f}}
    \end{align}
    Here, the introduction of dimensionless variables $\bar{\mathbf{u}}$, $\bar t$, and $\bar p$ has reduced the complexity of the original equations by consolidating physical parameters into dimensionless groups like the Reynolds Number (Re).

2.  **Identifying Dominant Physical Effects in Specific Regimes**

    Non-dimensionalization reveals the relative importance of various physical phenomena within a system.
    By examining the dimensionless numbers that emerge from the process, one can determine which terms in the equations are significant and which can be neglected under certain conditions.
    This identification is crucial for developing approximate solutions and understanding the behavior of the system in different regimes.

    For example, in high Reynolds number flows (Re $\gg 1$), inertial forces dominate over viscous forces, simplifying the Navier-Stokes equations by reducing the influence of the viscous term.
    Conversely, in low Reynolds number flows (Re $\ll 1$), viscous forces are predominant, and inertial terms can often be neglected.

+++

### Key Dimensionless Numbers

Several dimensionless numbers play pivotal roles in fluid dynamics and the study of PDEs.
Among the most important are the Reynolds Number (Re), Mach Number (Ma), and Prandtl Number (Pr).
Each of these numbers encapsulates the ratio of different physical effects, providing insight into the system's behavior.

+++

#### Reynolds Number (Re): Ratio of Inertial to Viscous Forces

The Reynolds Number is defined as:
\begin{align}
\text{Re} = \frac{\rho U L}{\mu} = \frac{U L}{\nu}
\end{align}
where:
* $\rho$ is the fluid density,
* $U$ is a characteristic velocity,
* $L$ is a characteristic length scale,
* $\mu$ is the dynamic viscosity,
* $\nu = \mu / \rho$ is the kinematic viscosity.

**Physical Interpretation:**

Reynolds Number quantifies the relative significance of inertial forces (associated with the fluid's motion) to viscous forces (associated with the fluid's internal friction). A high Re indicates that inertial forces dominate, leading to turbulent flow, while a low Re suggests that viscous forces are more influential, resulting in laminar flow.

**Applications:**

* **Flow Regimes:** Determining whether a flow will be laminar or turbulent based on Re.
* **Scale Modeling:** Ensuring dynamic similarity in wind tunnel experiments by matching Re between models and real-world scenarios.
* **Astrophysical Flows:** In stellar interiors, high Re can lead to turbulent convection, influencing energy transport and stellar evolution.

+++

#### Mach Number (Ma): Ratio of Flow Velocity to Speed of Sound

The Mach Number is defined as:
\begin{align}
\text{Ma} = \frac{U}{c}
\end{align}
where:
* $U$ is the characteristic flow velocity,
* $c$ is the speed of sound in the fluid.

**Physical Interpretation:**

Mach Number measures the compressibility effects in a flow. When Ma $< 1$, the flow is subsonic, and compressibility effects are negligible. When Ma $\approx 1$, the flow is transonic, and compressibility becomes significant. For Ma $> 1$, the flow is supersonic, and shock waves may form.

**Applications:**

* **Aerodynamics:** Designing aircraft and rockets by analyzing flow behavior at different Mach regimes.
* **Astrophysics:** Studying phenomena like shock waves in supernova explosions and supersonic jets from active galactic nuclei.
* **Explosive Events:** Understanding the propagation of shock waves generated by explosions.

+++

#### Prandtl Number (Pr): Ratio of Momentum Diffusivity to Thermal Diffusivity

The Prandtl Number is defined as:
\begin{align}
\text{Pr} = \frac{\nu}{\alpha} = \frac{\mu / \rho}{k / (\rho c_p)} = \frac{c_p \mu}{k}
\end{align}
where:
* $\nu$ is the kinematic viscosity,
* $\alpha = \frac{k}{\rho c_p}$ is the thermal diffusivity,
* $k$ is the thermal conductivity,
* $c_p$ is the specific heat at constant pressure.

**Physical Interpretation:**

Prandtl Number compares the rate at which momentum diffuses through a fluid to the rate at which heat diffuses.
A low Pr indicates that heat diffuses rapidly compared to momentum, while a high Pr suggests that momentum diffuses more quickly than heat.

**Applications:**

* **Heat Transfer:** Designing heat exchangers and cooling systems by understanding the relative rates of heat and momentum transfer.
* **Boundary Layer Analysis:** Determining the thickness of thermal and velocity boundary layers in fluid flows.
* **Astrophysical Processes:** Modeling heat conduction in stellar atmospheres where Pr influences energy transport mechanisms.

+++

### Applications of Dimensionless Numbers

Dimensionless numbers like Re, Ma, and Pr are instrumental in analyzing and predicting the behavior of physical systems across various domains:

1.  **Flow Control and Design**

    In engineering, these numbers guide the design of systems involving fluid flow.
    For example, ensuring that aircraft operate efficiently at desired Reynolds and Mach numbers is crucial for performance and safety.

2.  **Scale Modeling and Similitude**

  In experimental studies, ensuring that the dimensionless numbers of a model match those of the real system (dynamic similarity) allows for accurate predictions and validations of theoretical models through physical experiments.

3.  **Astrophysical Fluid Systems**

    In astrophysics, dimensionless numbers help in scaling and understanding fluid behaviors in vastly different environments:
    * **Stellar Convection:** High Reynolds numbers indicate turbulent convection currents within stars, affecting energy transport and mixing of stellar material.
    * **Accretion Disks:** Mach numbers determine the compressibility of gas in accretion disks around black holes, influencing the formation of shock waves and jet structures.
    * **Interstellar Medium:** Prandtl numbers aid in modeling the thermal and momentum diffusion in the diffuse interstellar gas, impacting star formation rates and the dynamics of molecular clouds.

+++

### Historical Anecdote: Fermi's Estimation of Atomic Bomb Energy Release

One of the most illustrative examples of the power of non-dimensionalization and dimensional analysis in physics is Enrico Fermi's remarkable estimation of the energy released during the first atomic bomb test, known as the Trinity test, conducted on July 16, 1945.
This anecdote not only highlights Fermi's ingenuity but also demostrate the practical utility of dimensionless numbers and scaling laws in tackling complex, real-world problems with limited empirical data.

During World War II, the Manhattan Project was a top-secret initiative aimed at developing nuclear weapons.
As the project progressed, scientists faced the daunting challenge of predicting the yield of the atomic bombs they were designing.
Precise calculations were hindered by the lack of comprehensive empirical data, necessitating innovative approaches to estimate explosive power accurately.

Enrico Fermi, a renowned physicist known for his ability to make quick, accurate estimates with minimal data, was tasked with predicting the energy release of the impending Trinity test.
His method exemplified dimensional analysis—a technique that uses the fundamental dimensions (such as mass, length, and time) to derive relationships between physical quantities.

Fermi's estimation process involved the following steps:
1.  **Observing the Shock Wave:**
   After the detonation, Fermi and his colleagues observed the speed at which the shock wave propagated through the air. Let's denote this speed as $U$.

2.  **Estimating the Shock Wave Radius:**
    By timing how long it took for the shock wave to reach a certain distance from the explosion site, Fermi estimated the radius $L$ of the blast wave after a specific duration $t$.

3.  **Applying the Sedov-Taylor Scaling:**
    The Sedov-Taylor solution describes the propagation of a strong shock wave from a point explosion in a homogeneous medium.
    According to this theory, the energy $E$ of the explosion is related to the shock wave's speed $U$, radius $L$, and the density $\rho$ of the surrounding air by the following relationship:
     \begin{align}
     E \sim \rho U^2 L^3
     \end{align}
     This equation is derived by balancing the kinetic energy imparted to the air by the shock wave with the energy of the explosion.

4.  **Calculating the Energy:**
    By substituting the estimated values of $\rho$, $U$, and $L$ into the above equation, Fermi arrived at an approximate value for the energy $E$ released during the explosion.

There is an arXiv paper on [this topic](https://arxiv.org/abs/2103.05784).
