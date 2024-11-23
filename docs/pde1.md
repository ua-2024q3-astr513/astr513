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

## Derivation of Fluid Dynamics Equations (25 minutes)

* Finite Volume Perspective
  * Overview of the finite volume method.
  * Derivation of fundamental fluid dynamics equations from conservation laws.

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
