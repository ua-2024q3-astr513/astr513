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

# Numerical Partial Differential Equation I: introduction

+++

## Introduction to Partial Differential Equations (PDEs)

Partial Differential Equations (PDEs) are fundamental tools in the mathematical modeling of various physical phenomena.
Unlike Ordinary Differential Equations (ODEs), which involve functions of a single variable and their derivatives, PDEs involve functions of multiple variables and their partial derivatives.
This distinction makes PDEs particularly powerful in describing systems where changes occur in more than one dimension, such as in space and time.

+++

### What are PDEs?

A Partial Differential Equation is an equation that relates the partial derivatives of a multivariable function.
In general form, a PDE can be written as:
\begin{align}
F\left(x_1, x_2, \ldots, x_n, u, \frac{\partial u}{\partial x_1}, \frac{\partial u}{\partial x_2}, \ldots, \frac{\partial^k u}{\partial x_1^{k_1} \partial x_2^{k_2} \ldots \partial x_n^{k_n}}\right) = 0
\end{align}
where $u = u(x_1, x_2, \ldots, x_n)$ is the unknown function, and $\partial u/\partial x_i$ denotes the partial derivatives of $u$ with respect to the variables $x_i$.

PDEs are essential in modeling continuous systems where the state of the system depends on multiple variables.
They appear in various fields such as physics, engineering, finance, and biology, describing phenomena like heat conduction, wave propagation, fluid dynamics, and quantum mechanics.

+++

### Definition and Significance in Modeling Continuous Systems

PDEs provide a framework for formulating problems involving functions of several variables and their rates of change.
They are indispensable in describing the behavior of physical systems where spatial and temporal variations are intrinsic.
For instance:

* **Heat Equation**: Models the distribution of heat (or temperature) in a given region over time.
  \begin{align}
  \frac{\partial u}{\partial t} = \alpha \nabla^2 u
  \end{align}

* **Wave Equation**: Describes the propagation of waves, such as sound or electromagnetic waves, through a medium.
  \begin{align}
  \frac{\partial^2 u}{\partial t^2} = c^2 \nabla^2 u
  \end{align}

- **Laplace's Equation**: Represents steady-state solutions where the system does not change over time, such as electric potential in a region devoid of charge.
  \begin{align}
  \nabla^2 u = 0
  \end{align}

The ability to model such diverse phenomena underscores the versatility and importance of PDEs in scientific and engineering disciplines.

+++

### Contrast with Ordinary Differential Equations (ODEs)

While both PDEs and ODEs involve differential operators, the key difference lies in the number of independent variables involved:

* **ODEs**: Involve functions of a single independent variable and their derivatives.
  They are typically used to model systems with dynamics that evolve over time without spatial variation.
  For example, the simple harmonic oscillator is described by the ODE:
  \begin{align}
  \frac{d^2 x}{dt^2} + \omega^2 x = 0
  \end{align}

* **PDEs**: Involve functions of multiple independent variables and their partial derivatives.
  They are used to model systems where changes occur across multiple dimensions, such as space and time.
  For example, the Navier-Stokes equations, which govern fluid flow, are a set of PDEs:
  \begin{align}
  \rho \left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} \right) = -\nabla p + \mu \nabla^2 \mathbf{u} + \mathbf{f}
  \end{align}
  where $\mathbf{u}$ is the velocity field, $p$ is the pressure, $\rho$ is the density, $\mu$ is the dynamic viscosity, and $\mathbf{f}$ represents body forces.

+++

### Application in Astrophysical Fluid Systems

PDEs play a pivotal role in modeling complex astrophysical phenomena, particularly those involving fluid dynamics.
Astrophysical fluid systems, such as stellar interiors, accretion disks around black holes, and interstellar gas clouds, exhibit intricate behaviors governed by the laws of fluid mechanics and thermodynamics. 

For example, the **Euler Equations** for inviscid flow and the **Magnetohydrodynamic (MHD) Equations** for conducting fluids are essential in understanding phenomena like solar flares and the dynamics of the interstellar medium.
These equations are inherently PDEs due to the multi-dimensional and time-dependent nature of the systems they describe.

Consider the modeling of **stellar convection**, where hot plasma rises and cooler plasma sinks within a star's interior.
This process is governed by the Navier-Stokes equations coupled with the heat equation, forming a system of PDEs that describe the velocity field, temperature distribution, and pressure variations within the star.

Furthermore, **gravitational waves**, ripples in spacetime caused by massive astrophysical events like neutron star mergers, are described by the Einstein Field Equations, which are a set of nonlinear PDEs.
Solving these equations numerically requires sophisticated techniques due to their complexity and the need for high precision in predictions.

In summary, PDEs are indispensable in astrophysics for modeling and understanding the dynamic and multifaceted nature of celestial phenomena.
The ability to solve these equations, either analytically or numerically, provides deep insights into the behavior and evolution of the universe.

+++

## Derivation of Fluid Dynamics Equations

Understanding the fundamental equations that govern fluid motion is essential for both theoretical studies and practical applications in engineering and astrophysics.
This section delves into the derivation of the primary fluid dynamics equations—namely, the Continuity, Momentum (Navier-Stokes), and Energy equations—using the integral forms of conservation laws and applying Green's and Stokes' theorems.
Additionally, we explore an alternative derivation from the Boltzmann Equation using the Moment Method, providing a comprehensive foundation for modeling fluid behavior in various contexts, including complex astrophysical systems.

+++

### Integral Forms of Conservation Laws

Fluid dynamics is fundamentally rooted in the principles of conservation of mass, momentum, and energy.
These principles can be expressed in their integral forms, which consider the behavior of fluid quantities within a finite control volume.

+++

1. **Conservation of Mass (Continuity Equation):**

   The integral form of the continuity equation states that the rate of change of mass within a control volume $V$ is equal to the net mass flux entering or leaving the volume through its boundary $\partial V$:
   \begin{align}
   \frac{d}{dt} \int_V \rho \, dV + \oint_{\partial V} \rho \mathbf{u} \cdot \mathbf{n} \, dS = 0
   \end{align}
   where:
   * $\rho$ is the fluid density,
   * $\mathbf{u}$ is the velocity field,
   * $\mathbf{n}$ is the outward-pointing unit normal vector on $\partial V$.

+++

2. **Conservation of Momentum (Navier-Stokes Equations):**

   The integral form of the momentum conservation law accounts for the forces acting on the fluid within the control volume.
   It can be expressed as:
   \begin{align}
   \frac{d}{dt} \int_V \rho \mathbf{u} \, dV + \oint_{\partial V} \rho \mathbf{u} (\mathbf{u} \cdot \mathbf{n}) \, dS = \oint_{\partial V} \mathbf{\Pi} \cdot \mathbf{n} \, dS + \int_V \rho \mathbf{f} \, dV
   \end{align}
   where:
   * $\mathbf{\Pi}$ is the stress tensor,
   * $\mathbf{f}$ represents body forces (e.g., gravity).

+++

3. **Conservation of Energy:**

   The integral form of the energy conservation law relates the rate of change of energy within the control volume to the net energy flux and work done by forces:
   \begin{align}
   \frac{d}{dt} \int_V \rho e \, dV + \oint_{\partial V} \rho e \mathbf{u} \cdot \mathbf{n} \, dS = \oint_{\partial V} \mathbf{q} \cdot \mathbf{n} \, dS + \oint_{\partial V} \mathbf{\Pi} \cdot \mathbf{u} \, dS + \int_V \rho \mathbf{f} \cdot \mathbf{u} \, dV
   \end{align}
   where:
   - $e$ is the specific internal energy,
   - $\mathbf{q}$ is the heat flux vector.

+++

### From Integral to Differential Forms Using Green's and Stokes' Theorems

To transition from the integral to the differential forms of these conservation laws, we employ Green's and Stokes' theorems, which relate volume integrals to surface integrals.

The Divergence Theorem converts a volume integral of a divergence of a vector field into a surface integral over the boundary of the volume:
\begin{align}
\oint_{\partial V} \mathbf{F} \cdot \mathbf{n} \, dS = \int_V \nabla \cdot \mathbf{F} \, dV
\end{align}

+++

#### Application to Conservation of Mass

Applying the Divergence Theorem to the continuity equation:
\begin{align}
\frac{d}{dt} \int_V \rho \, dV + \oint_{\partial V} \rho \mathbf{u} \cdot \mathbf{n} \, dS = 0 \quad 
\Rightarrow \quad \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0
\end{align}

This yields the **Continuity Equation** in differential form:
\begin{align}
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0
\end{align}

+++

#### Application to Conservation of Momentum

Applying the Divergence Theorem to the momentum equation:
\begin{align}
\frac{d}{dt} \int_V \rho \mathbf{u} \, dV + \oint_{\partial V} \rho \mathbf{u} (\mathbf{u} \cdot \mathbf{n}) \, dS
&= \oint_{\partial V} \mathbf{\Pi} \cdot \mathbf{n} \, dS + \int_V \rho \mathbf{f} \, dV\\
\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u})
&= \nabla \cdot \mathbf{\Pi} + \rho \mathbf{f}
\end{align}

To obtain the **Navier-Stokes Equations**, we express the stress tensor $\mathbf{\Pi}$ for a Newtonian fluid:
\begin{align}
\mathbf{\Pi} = -p \mathbf{I} + \mu \left( \nabla \mathbf{u} + (\nabla \mathbf{u})^T \right) + \lambda (\nabla \cdot \mathbf{u}) \mathbf{I}
\end{align}
where:
* $p$ is the pressure,
* $\mu$ is the dynamic viscosity,
* $\lambda$ is the second coefficient of viscosity,
* $\mathbf{I}$ is the identity tensor.

Substituting $\mathbf{\Pi}$ into the momentum equation and simplifying under the assumption of incompressible flow ($ \nabla \cdot \mathbf{u} = 0$) leads to:
\begin{align}
\rho \left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} \right) = -\nabla p + \mu \nabla^2 \mathbf{u} + \rho \mathbf{f}
\end{align}

+++

#### Application to Conservation of Energy

Similarly, applying the Divergence Theorem to the energy equation:
\begin{align}
\frac{d}{dt} \int_V \rho e \, dV + \oint_{\partial V} \rho e \mathbf{u} \cdot \mathbf{n} \, dS 
&= \oint_{\partial V} \mathbf{q} \cdot \mathbf{n} \, dS + \oint_{\partial V} \mathbf{\Pi} \cdot \mathbf{u} \, dS + \int_V \rho \mathbf{f} \cdot \mathbf{u} \, dV \\
\frac{\partial (\rho e)}{\partial t} + \nabla \cdot (\rho e \mathbf{u}) 
&= \nabla \cdot \mathbf{q} + \mathbf{\Pi} : \nabla \mathbf{u} + \rho \mathbf{f} \cdot \mathbf{u}
\end{align}

Assuming Fourier's law for heat conduction ($\mathbf{q} = -k \nabla T$) and substituting the expression for $ \mathbf{\Pi}$, we obtain the **Energy Equation** in differential form:
\begin{align}
\frac{\partial (\rho e)}{\partial t} + \nabla \cdot (\rho e \mathbf{u}) = \nabla \cdot (k \nabla T) + \Phi + \rho \mathbf{f} \cdot \mathbf{u}
\end{align}
where $\Phi$ represents the viscous dissipation function.

+++

* From Boltzmann to Navier-Stokes
  * Boltzmann Equation Overview
    * Introduction to the particle distribution function.
  * Moment Method:
    * Derivation of the Continuity Equation.
    * Derivation of the Momentum Equation.
    * Derivation of the Energy Equation.
  * Assumptions & Physical Meaning
  * Viscosity, pressure gradients, and other key terms.
  * Limiting cases and their physical interpretations.

* Significance in Fluid Dynamics

* Role of PDEs in governing real-world flows (e.g., airflow, ocean currents).

+++

## Classification of Partial Differential Equations (PDEs) (10 minutes)

* Types of PDEs:
  * Elliptic: Example – Laplace's equation.
  * Parabolic: Example – Heat equation.
  * Hyperbolic: Example – Wave equation.

* Physical Intuition

* Understanding the nature and applications of each type.

* Examples of physical phenomena governed by each class of PDEs.

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

## Numerical Techniques for Solving PDEs (25 minutes)

* Finite Difference Methods (FDM)
  * Forward Time Centered Space (FTCS):
    * Formulation and inherent instability issues.
  * Von Neumann Stability Analysis
    * Defining the amplification factor.
    * Stability criteria.
  * Stabilized Methods
    * Lax Method:
      * Introduction of numerical dissipation for stability.
    * Courant-Friedrichs-Lewy (CFL) Condition:
      * Ensuring numerical stability through timestep restrictions.

* Method Comparison
  * Evaluating methods based on stability, accuracy, and computational cost.
  * Practical considerations for method selection in various scenarios.

+++

## Practical Example: The Heat Equation (5 minutes)

* 1D Heat Equation Overview
  * Derivation from conservation principles: $u_t = \alpha u_{xx}$.
  * Analytical solution using separation of variables.
  * Discussion on boundary and initial conditions.

* Numerical Solution Approach
  * Basic finite difference implementation steps.
  * Impact and importance of adhering to the CFL condition for stability.
