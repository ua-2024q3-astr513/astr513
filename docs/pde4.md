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

# Numerical Partial Differential Equation IV: Spectral Methods

+++

## **Introduction to Spectral Methods**

Spectral methods are a powerful class of numerical techniques used to solve partial differential equations (PDEs) by representing the solution as a global expansion of orthogonal basis functions.
Unlike finite difference or finite volume methods, which approximate the solution locally at discrete grid points, spectral methods leverage the smoothness of the solution to achieve high accuracy with fewer degrees of freedom.
This makes spectral methods particularly suitable for problems involving smooth and periodic solutions, such as fluid dynamics and wave propagation.

+++

**Why Use Spectral Methods?**

1. **High Accuracy:**
   Spectral methods are renowned for their ability to achieve exponential convergence rates for smooth solutions.
   By representing the solution as a sum of global basis functions, spectral methods resolve fine details with far fewer grid points compared to traditional numerical methods.

2. **Efficient Computation:**
   Operations such as differentiation and integration become algebraic manipulations in spectral space.
   This allows for fast and efficient computation, especially when coupled with Fast Fourier Transforms (FFT).

3. **Energy-Conserving Properties:**
   By carefully choosing the basis functions and truncating higher-order modes (Galerkin truncation), spectral methods naturally conserve key physical quantities like energy and enstrophy.
   This property is critical in applications such as turbulence modeling and geophysical fluid dynamics.

4. **Applications in Science and Engineering:**
   Spectral methods are widely used in areas such as:
   * Astrophysical fluid dynamics.
   * Weather prediction and ocean modeling.
   * Quantum mechanics and wave phenomena.

+++

**Core Idea of Spectral Methods**

The central idea of spectral methods is to approximate a function $u(x)$ as a series of orthogonal basis functions.
For example, in the Fourier spectral method, the function is expanded as:
\begin{align}
u(x) = \sum_{k=-N/2}^{N/2} \hat{u}_k e^{i k x},
\end{align}
where $\hat{u}_k$ are the Fourier coefficients representing the contribution of each basis function.
The advantage of this representation is that operations like differentiation transform into simple multiplications in spectral space.
For example:
\begin{align}
\frac{\partial u}{\partial x} = \sum_{k=-N/2}^{N/2} i k \hat{u}_k e^{i k x}.
\end{align}
This property drastically simplifies the computation of derivatives, making spectral methods particularly attractive for solving PDEs.

+++

**Periodic vs. Non-Periodic Domains**

Spectral methods are naturally suited to periodic domains, where Fourier basis functions $e^{i k x}$ form an orthogonal basis.
For non-periodic problems, alternative basis functions like Chebyshev polynomials or Legendre polynomials are used.
These basis functions maintain the accuracy of spectral methods while adapting to non-periodic boundary conditions.

+++

**Strengths**
* **Exponential Convergence:** For smooth problems, spectral methods outperform traditional numerical methods in terms of accuracy.
* **Global Representation:** Captures global features of the solution with fewer degrees of freedom.
* **Efficiency with FFT:** Fast Fourier Transforms enable rapid computation of spectral coefficients.

**Limitations**
* **Smoothness Requirement:** Spectral methods perform poorly for non-smooth or discontinuous solutions, where the Gibbs phenomenon introduces oscillations.
* **Complexity for Non-Periodic Domains:** Implementing spectral methods for non-periodic problems requires specialized basis functions and quadrature rules.
* **Global Coupling:** Each mode influences the entire domain, which can lead to computational challenges for very large systems.

+++

**Spectral Methods in Fluid Dynamics**

In fluid dynamics, spectral methods are commonly used to solve the incompressible Navier-Stokes equations.
By transforming the equations into spectral space, the nonlinear terms can be computed efficiently while maintaining conservation properties.
For incompressible flows, the vorticity-streamfunction formulation is particularly suitable, as it eliminates the pressure term and reduces the computational complexity.

In this lecture, we will focus on applying spectral methods to solve the **2D incompressible hydrodynamics equations**. Specifically, we will:
1. Derive the governing equations in spectral space.
2. Discuss the vorticity-streamfunction formulation and its advantages.
3. Implement a spectral solver to simulate the evolution of vorticity in a 2D periodic domain.

This introduction sets the stage for understanding how spectral methods leverage mathematical elegance and computational efficiency to solve complex PDEs with remarkable accuracy.

+++

## Introduction to 2D Incompressible Hydrodynamics

The hydrodynamic equations govern the conservation of mass and momentum in a fluid.
In their compressible form (that we derived in previous lectures), they are written as:
\begin{align}
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) &= 0, \quad \text{(Continuity Equation)} \\
\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \mathbf{u}) &= -\nabla p + \mu \nabla^2 \mathbf{u} + \mathbf{f}, \quad \text{(Momentum Equation)}
\end{align}
where:
* $\rho$ is the density,
* $\mathbf{u}$ is the velocity field,
* $p$ is the pressure,
* $\mu$ is the dynamic viscosity,
* $\mathbf{f}$ is an external force.

+++

In the **incompressible limit**, the sound speed approaches infinite $c \rightarrow \infty$.
For simplicity, the density $\rho$ can be assumed constant, and the continuity equation reduces to the **incompressibility condition**:
\begin{align}
\nabla \cdot \mathbf{u} = 0.
\end{align}
Substituting this condition into the momentum equation simplifies the Navier-Stokes equations to:
\begin{align}
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = -\nabla p + \nu \nabla^2 \mathbf{u} + \mathbf{f},
\end{align}
where $\nu = \mu / \rho$ is the kinematic viscosity.
These equations describe the flow of incompressible fluids and are widely used in modeling small-scale laboratory experiments and large-scale geophysical flows.

+++

While the incompressible Navier-Stokes equations also apply to three-dimensional flows, many physical systems can be effectively approximated as two-dimensional.
For example:
* Atmospheric flows and ocean currents are largely horizontal due to their vast spatial extent compared to their depth.
* Thin liquid films and confined flows are geometrically restricted to two dimensions.

In 2D, the dynamics exhibit unique features that distinguish them from 3D flows:

**Conservation of Enstrophy**

In 2D, the vorticity $w = \nabla \times \mathbf{u}$ is a scalar field.
Its evolution is governed by the **vorticity transport equation**, which conserves both **energy** and **enstrophy** in the absence of dissipation:
\begin{align}
E &= \frac{1}{2} \int |\mathbf{u}|^2 \, dx \, dy\\
Z &= \frac{1}{2} \int w^2 \, dx \, dy.
\end{align}

Energy conservation governs the total kinetic energy of the system, while enstrophy conservation introduces a second constraint that strongly influences the flow dynamics.

+++

**Inverse Energy Cascade**

A striking feature of 2D turbulence is the **inverse energy cascade**.
In 3D turbulence, energy flows from large scales (low wavenumbers) to small scales (high wavenumbers) and is dissipated by viscosity.
In 2D, however, energy flows in the opposite direction, from small scales to large scales, leading to the formation of large, coherent structures like cyclones and anticyclones.
This behavior is directly tied to the dual conservation of energy and enstrophy.

+++

## Simplified Governing Equations for 2D Flows

To simplify the analysis of incompressible flows, we introduce the **streamfunction** $\psi$, which ensures the incompressibility condition is automatically satisfied:
\begin{align}
\mathbf{u} = \nabla \times (\psi \mathbf{\hat{z}}), \quad u_x = \frac{\partial \psi}{\partial y}, \quad u_y = -\frac{\partial \psi}{\partial x}.
\end{align}

The vorticity $w$, defined as $w = \nabla \times \mathbf{u}$, relates to $\psi$ via the **Poisson equation**:
\begin{align}
w = -\nabla^2 \psi.
\end{align}

The vorticity transport equation describes the evolution of $w$:
\begin{align}
\frac{\partial w}{\partial t} + \mathbf{u} \cdot \nabla w = \nu \nabla^2 w + \mathbf{f}_w,
\end{align}

where $\mathbf{f}_w$ is a forcing term that may drive the flow. This formulation eliminates the pressure term and reduces the computational complexity, making it ideal for numerical simulations using spectral methods.
