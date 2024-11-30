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

# Numerical Partial Differential Equation III: Finite Volume

+++

In the previous lecture, we explored three fundamental finite difference schemes used to solve the linear advection equation: the Forward Time Centered Space (FTCS) scheme, the Upwind scheme, and the Lax-Wendroff scheme.
Each of these methods offers a unique balance between simplicity, accuracy, and stability.
The FTCS scheme, while straightforward to implement, was found to be inherently unstable for advection problems.
The Upwind scheme addressed some stability issues by introducing numerical diffusion, which helped in suppressing non-physical oscillations but at the cost of smearing sharp gradients.
The Lax-Wendroff scheme provided a higher order of accuracy by incorporating both first and second-order spatial derivatives, reducing numerical diffusion and better capturing wave propagation, though it introduced dispersive errors.

Despite their effectiveness in one-dimensional problems, finite difference methods (FDM) like FTCS, Upwind, and Lax-Wendroff encounter significant challenges when applied to more complex scenarios.
Specifically, FDMs struggle with handling complex geometries and discontinuities, such as shock waves, which are common in multi-dimensional fluid dynamics problems.
Additionally, FDMs often face difficulties in conserving physical quantities, such as mass, momentum, and energy, across control volumes, especially when dealing with irregular meshes or boundaries.

To overcome these limitations, **Finite Volume Methods (FVM)** have been developed as a robust alternative to finite difference approaches.
Unlike FDMs, which focus on approximating derivatives at discrete points, FVMs emphasize the conservation of physical quantities within discrete control volumes.
This integral approach ensures that fluxes entering and leaving a control volume are balanced, inherently satisfying conservation laws.
Consequently, FVMs are well-suited for solving complex, multi-dimensional problems involving sharp gradients and discontinuities, making them a preferred choice in computational fluid dynamics and related fields.

In this lecture, we will delve into the motivation behind Finite Volume Methods, understand their foundational principles, and explore how traditional finite difference schemes can be reformulated within the FVM framework.
We will also introduce a simple finite volume method tailored for solving two-dimensional shock problems, highlighting the advantages of FVMs in handling complex flow phenomena.

+++

## Motivation for Finite Volume Methods

Finite Volume Methods (FVM) have become a cornerstone in computational fluid dynamics and other areas involving the simulation of physical phenomena governed by conservation laws.
The motivation behind adopting FVM over traditional Finite Difference Methods (FDM) stems from several key advantages that address the limitations inherent in FDMs.

1. **Conservation Laws and Fluxes**

   At the heart of many physical systems are conservation laws, which dictate that certain quantities such as mass, momentum, and energy remain constant within a closed system.
   FVM inherently respects these conservation principles by integrating the governing equations over discrete control volumes.
   Unlike FDMs, which approximate derivatives at specific points, FVMs focus on the fluxes of conserved quantities across the boundaries of control volumes.
   By ensuring that the net flux into a control volume equals the rate of change of the conserved quantity within that volume, FVMs maintain local conservation properties.
   This integral approach guarantees that the discretized equations faithfully represent the continuous conservation laws, leading to more accurate and physically consistent solutions.

3. **Handling Discontinuities and Shocks**

   One of the significant challenges in numerical simulations is accurately capturing discontinuities and shock waves, which are abrupt changes in the flow properties.
   FDMs, especially lower-order schemes like FTCS, often struggle with these features due to excessive numerical diffusion, which smears out sharp gradients and can distort the physical phenomena being modeled.
   In contrast, FVMs are adept at handling such discontinuities.
   The flux-based formulation allows FVMs to more precisely track the movement and interaction of shocks and contact discontinuities.
   By accurately computing the fluxes at control volume interfaces, FVMs can maintain the integrity of sharp gradients without introducing significant numerical artifacts.
   This capability makes FVMs particularly suitable for simulating high-speed flows and other scenarios where shock waves play a crucial role.

5. **Advantages of FVM Over FDM and Finite Element Methods (FEM)**

   Finite Volume Methods offer several advantages over both Finite Difference Methods and Finite Element Methods:

   * **Local Conservation Properties:** FVMs ensure that the conserved quantities are maintained within each control volume.
     This local conservation is essential for accurately modeling physical systems, especially over long simulation times or in complex flow conditions.

   * **Flexibility in Handling Complex Geometries and Boundary Conditions:** FVMs are inherently more flexible when dealing with irregular geometries and complex boundary conditions.
     The ability to define control volumes of various shapes and sizes allows FVMs to conform to intricate domain boundaries, making them suitable for real-world applications where geometrical complexity is common.

   * **Compatibility with Unstructured Meshes:** FVMs seamlessly integrate with unstructured meshes, which are essential for discretizing domains with irregular shapes.
     Unstructured meshes allow for greater adaptability and refinement in regions requiring higher accuracy, such as areas with steep gradients or complex flow features.
     This compatibility enhances the overall accuracy and efficiency of the numerical simulations.

   * **Ease of Incorporating Physical Models:** FVMs facilitate the incorporation of additional physical models, such as turbulence models or multi-phase flow dynamics, within the control volume framework.
     This modularity enables the extension of FVMs to a wide range of applications beyond simple advection or diffusion processes.

   These advantages make FVMs a powerful tool for accurately and efficiently solving a broad spectrum of physical problems that are challenging for traditional Finite Difference Methods and Finite Element Methods.

+++

## Basics of Finite Volume Methods

Finite Volume Methods (FVM) are widely used for solving Partial Differential Equations (PDEs) that express conservation laws.
Unlike Finite Difference Methods (FDM), which approximate derivatives at discrete points, FVM focuses on conserving quantities within discrete control volumes.
This section introduces the foundational concepts of FVM, including control volumes, integral formulation of conservation laws, discretization of fluxes, and numerical flux functions.

+++

### Control Volumes and Integral Formulation

At the core of Finite Volume Methods is the concept of dividing the computational domain into small, non-overlapping regions called **control volumes**.
Each control volume encompasses a portion of the physical domain and is bounded by its interfaces with neighboring control volumes.

**Definition of Control Volumes:**
* A **control volume** is a finite region in space over which conservation laws are applied.
* The entire domain is partitioned into these control volumes, ensuring that they cover the domain without gaps or overlaps.

+++

**Integral Form of Conservation Laws:**
Finite Volume Methods are based on the integral form of conservation laws. Consider a general conservation equation:
\begin{align}
\frac{\partial \phi}{\partial t} + \nabla \cdot \mathbf{F} = S
\end{align}
where:
- $\phi$ is the conserved quantity (e.g., mass, momentum, energy),
- $\mathbf{F}$ is the flux vector representing the flow of $\phi$,
- $S$ is a source term.

+++

**Integral Formulation:**
Integrate the conservation equation over a control volume $V$:
\begin{align}
\int_V \frac{\partial \phi}{\partial t} \, dV + \int_{V} \nabla\cdot\mathbf{F} \, dV = \int_V S \, dV
\end{align}

+++

**Applying the Divergence Theorem:**
The integral of the divergence of $\mathbf{F}$ over the control volume is converted to a surface integral:
\begin{align}
\frac{d}{dt} \int_V \phi \, dV + \int_{\partial V} \mathbf{F} \cdot \mathbf{n} \, dS = \int_V S \, dV
\end{align}
where $\mathbf{n}$ is the outward-pointing unit normal vector on the boundary $\partial V$ of volume $V$.
This equation states that the rate of change of $\phi$ within the control volume plus the net flux of $\phi$ across its boundaries equals the total source of $\phi$ within the volume.

Note that this is the same formulation that we derived the fluid dynamic equations using conservation laws.

+++

### Discretization of Fluxes

To solve the integral conservation laws numerically, the flux integrals across the control volume boundaries must be approximated.
The discretization process involves several steps:

1. **Mesh Generation:**
   * Divide the computational domain into a finite number of control volumes.
   * Common mesh types include structured (e.g., rectangular, triangular) and unstructured meshes.
     
2. **Control Volume Centroid:**
   * Identify the centroid or a representative point within each control volume, typically denoted as $(x_i, y_i)$ in 2D or $(x_i, y_i, z_i)$ in 3D.

3. **Flux Approximation:**
   * Approximate the flux $\mathbf{F}$ at each face of the control volume.
   * The flux at a face is determined by evaluating the flux function based on the values of $\phi$ from adjacent control volumes.

4. **Surface Integral Discretization:**
   * Replace the continuous surface integral with a discrete sum over all faces of the control volume:
     \begin{align}
     \int_{\partial V} \mathbf{F} \cdot \mathbf{n} \, dS \approx \sum_{\text{faces}} \mathbf{F}_{\text{face}} \cdot A_{\text{face}}
     \end{align}
     where $A_{\text{face}}$ is the length (for 2D problem) or area (for 3D problem) of the "face".

+++

**Example for One-Dimensional Advection:**
Consider the linear advection equation in one dimension:
\begin{align}
\frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0
\end{align}
For a control volume $[x_{i-\frac{1}{2}}, x_{i+\frac{1}{2}}]$, the integral form becomes:
\begin{align}
\frac{d}{dt} \int_{x_{i-\frac{1}{2}}}^{x_{i+\frac{1}{2}}} u \, dx + \left[ c u \right]_{x_{i-\frac{1}{2}}}^{x_{i+\frac{1}{2}}} = 0
\end{align}
Discretizing the fluxes at the boundaries:
\begin{align}
\frac{d u_i}{dt} \Delta x + c (u_{i+\frac{1}{2}} - u_{i-\frac{1}{2}}) = 0
\end{align}
Solving for $\frac{d u_i}{dt}$:
\begin{align}
\frac{d u_i}{dt} = -\frac{c}{\Delta x} (u_{i+\frac{1}{2}} - u_{i-\frac{1}{2}})
\end{align}

+++

### Numerical Flux Functions

The accuracy and stability of the Finite Volume Method heavily depend on the choice of **numerical flux functions**.
Numerical flux functions approximate the physical fluxes at the interfaces between control volumes based on the values of $\phi$ from adjacent cells.
Several types of numerical flux functions exist, each with its own advantages and applicability:

1. **Central Flux:**
   * A simple average of the fluxes from adjacent control volumes.
   * **Formula:**
     \begin{align}
     \mathbf{F}_{i+\frac{1}{2}} = \frac{1}{2} \left( \mathbf{F}(u_i) + \mathbf{F}(u_{i+1}) \right)
     \end{align}
   * **Advantages:**
     * Easy to implement.
     * Second-order accurate for smooth solutions.
   * **Disadvantages:**
     * Can be unstable for hyperbolic equations like advection.
     * Does not account for wave direction, leading to numerical oscillations.

2. **Upwind Flux:**
   * Takes the flux from the upstream direction based on the wave speed $c$.
   * **Formula for One-Dimensional Advection:**
     \begin{align}
     \mathbf{F}_{i+\frac{1}{2}} = 
     \begin{cases}
       \mathbf{F}(u_i) & \text{if } c > 0 \\
       \mathbf{F}(u_{i+1}) & \text{if } c < 0
     \end{cases}
     \end{align}
   * **Advantages:**
     * Enhances stability by considering wave direction.
     * Reduces numerical oscillations near discontinuities.
   * **Disadvantages:**
     * Introduces numerical diffusion, which can smear sharp gradients.

3. **Riemann Solver-Based Fluxes:**
   * Solve the Riemann problem at each interface to determine the flux.
   * **Examples:**
     * Exact Riemann Solver.
     * Roe's Approximate Riemann Solver.
   * **Advantages:**
     * Accurately captures wave interactions.
     * Suitable for complex systems of conservation laws.
   * **Disadvantages:**
     * More computationally intensive.
     * Requires solving characteristic equations.

4. **HLL and HLLC Fluxes:**
   * Variants of Riemann solvers that simplify the wave structure.
   * **HLL (Harten, Lax, and van Leer):**
     * Considers only the fastest left and right waves.
   * **HLLC (HLL with Contact):**
     * Includes contact discontinuities in addition to shock waves.
   * **Advantages:**
     * Balances accuracy and computational efficiency.
     * Handles a wide range of flow conditions.
   * **Disadvantages:**
     * Still more complex than central or upwind fluxes.

**Choosing the Appropriate Numerical Flux Function:**
The selection of a numerical flux function depends on the specific requirements of the problem:
- **For Simple Problems:** Central or upwind fluxes may suffice.
- **For Problems with Shocks and Complex Wave Interactions:** Riemann solver-based fluxes or HLL/HLLC fluxes are more appropriate.

**Example: Upwind Flux in One-Dimensional Advection**
For the linear advection equation with $c > 0$, the upwind flux at interface $i+\frac{1}{2}$ is:
\begin{align}
F_{i+\frac{1}{2}} = c u_i
\end{align}
Substituting into the finite volume discretization:
\begin{align}
\frac{d u_i}{dt} = -\frac{c}{\Delta x} (u_i - u_{i-1})
\end{align}
This corresponds to the Upwind finite difference scheme when formulated within the finite volume framework.

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

# Physical Parameters
gamma = 1.4  # Specific heat ratio for air

# Computational Domain
L  = 1.0    # Length of the domain (meters)
N  = 100    # Number of cells
dx = L / N  # Cell width
x  = np.linspace(0.5*dx, L - 0.5*dx, N)  # Cell centers

# Time-Stepping Parameters
CFL     = 0.5   # Courant number
t_final = 0.25  # Final time (seconds)
t       = 0.0   # Initial time
```

```{code-cell} ipython3
# Initialize Conserved Variables [rho, rho*u, E]
U = np.zeros((N, 3))

# Initial Conditions
rho0 = np.where(x < 0.5, 1.0, 0.125)          # Density
u0   = np.zeros(N)                            # Velocity
p0   = np.where(x < 0.5, 1.0, 0.1)            # Pressure
E0   = p0 / (gamma - 1) + 0.5 * rho0 * u0**2  # Total Energy

U[:, 0] = rho0       # Density
U[:, 1] = rho0 * u0  # Momentum
U[:, 2] = E0         # Energy
```

```{code-cell} ipython3
# Function to convert conserved to primitive variables
def conserved_to_primitive(U):
    rho = U[...,0]
    u   = U[...,1] / rho
    E   = U[...,2]
    p   = (gamma - 1) * (E - 0.5 * rho * u**2)
    return rho, u, p

# HLL Riemann Solver
def HLL_flux(UL, UR):

    rhoL, uL, pL = conserved_to_primitive(UL)
    rhoR, uR, pR = conserved_to_primitive(UR)
    EL = UL[2]
    ER = UR[2]
    cL = np.sqrt(gamma * pL / rhoL)
    cR = np.sqrt(gamma * pR / rhoR)
    
    # Wave Speeds
    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)
    
    # Fluxes for left and right states
    FL = np.array([
        rhoL * uL,
        rhoL * uL**2 + pL,
        uL   * (EL + pL)
    ])
    
    FR = np.array([
        rhoR * uR,
        rhoR * uR**2 + pR,
        uR   * (ER + pR)
    ])
    
    # HLL Flux Calculation
    if SL > 0:
        return FL
    elif SL <= 0 and SR >= 0:
        return (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL)
    else:
        return FR
```

```{code-cell} ipython3

```

```{code-cell} ipython3
# Time-Stepping Loop
time_steps = [t]
U_history  = [U.copy()]

while t < t_final:
    rho, u, p = conserved_to_primitive(U)
    c     = np.sqrt(gamma * p / rho)
    s_max = np.max(np.abs(u) + c)
    dt    = CFL * dx / s_max
    if t + dt > t_final:
        dt = t_final - t
    
    # Compute Fluxes at Interfaces
    flux = np.zeros((N+1, 3))
    for i in range(1,N):
        flux[i] = HLL_flux(U[i-1], U[i])

    # Boundary conditions on flux
    flux[ 0] = flux[ 1]
    flux[-1] = flux[-2]
    
    # Update Conserved Variables
    U -= (dt / dx) * (flux[1:] - flux[:-1])

    # Update Time and Variables
    t += dt
    time_steps.append(t)
    U_history.append(U.copy())
```

```{code-cell} ipython3
# Visualization of Results
rho, u, p = conserved_to_primitive(U)

plt.figure(figsize=(14, 4))

plt.subplot(1, 3, 1)
plt.plot(x, rho0, ':', label='Density', color='C0')
plt.plot(x, rho,       label='Density', color='C0')
plt.xlabel('Position')
plt.ylabel('Density')
plt.title('Density at t = {:.3f} s'.format(t_final))
plt.legend()
plt.grid()

plt.subplot(1, 3, 2)
plt.plot(x, u0, ':', label='Velocity', color='C1')
plt.plot(x, u,       label='Velocity', color='C1')
plt.xlabel('Position')
plt.ylabel('Velocity')
plt.title('Velocity at t = {:.3f} s'.format(t_final))
plt.legend()
plt.grid()

plt.subplot(1, 3, 3)
plt.plot(x, p0, ':', label='Pressure', color='C2')
plt.plot(x, p,       label='Pressure', color='C2')
plt.xlabel('Position')
plt.ylabel('Pressure')
plt.title('Pressure at t = {:.3f} s'.format(t_final))
plt.legend()
plt.grid()

plt.tight_layout()
```