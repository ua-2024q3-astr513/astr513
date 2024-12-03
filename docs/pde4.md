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
