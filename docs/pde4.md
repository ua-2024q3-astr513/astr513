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

## Vorticity-Streamfunction Formulation

To simplify the mathematical and computational treatment of 2D incompressible hydrodynamics, the governing equations are often reformulated in terms of the **vorticity** $w$ and the **streamfunction** $\psi$.
This formulation has several advantages: it eliminates the pressure term from the equations, reduces the number of variables, and ensures incompressibility is automatically satisfied.
In this section, we derive the vorticity-streamfunction formulation, define the Jacobian determinant to handle the nonlinear advection term, and introduce additional physical effects such as Ekman damping and the beta-plane approximation.

+++

### Definitions and Key Relationships

For 2D incompressible flows, the velocity field $\mathbf{u} = (u_x, u_y)$ can be expressed in terms of a scalar function, the **streamfunction** $\psi(x, y, t)$, as:
\begin{align}
\mathbf{u} = \nabla \times (\psi \mathbf{\hat{z}}),
\end{align}
where $\mathbf{\hat{z}}$ is the unit vector perpendicular to the 2D plane.
In component form:
\begin{align}
u_x = \frac{\partial \psi}{\partial y}, \quad u_y = -\frac{\partial \psi}{\partial x}.
\end{align}

This representation automatically satisfies the incompressibility condition:
\begin{align}
\nabla \cdot \mathbf{u} = \frac{\partial u_x}{\partial x} + \frac{\partial u_y}{\partial y} = 0.
\end{align}

+++

The vorticity $w$ is a scalar quantity in 2D, defined as the curl of the velocity field:
\begin{align}
w = \nabla \times \mathbf{u}.
\end{align}

Using the velocity components in terms of the streamfunction, the vorticity can be written as:
\begin{align}
w = \frac{\partial u_y}{\partial x} - \frac{\partial u_x}{\partial y} = -\nabla^2 \psi.
\end{align}

Thus, the vorticity and streamfunction are related by the **Poisson equation**:
\begin{align}
w = -\nabla^2 \psi.
\end{align}

This relationship allows the velocity field to be computed from the vorticity by solving the Poisson equation for $\psi$, and then taking derivatives of $\psi$ to find $u_x$ and $u_y$.

+++

### Governing Equation for Vorticity

The vorticity transport equation is derived from the incompressible Navier-Stokes equations.
Taking the curl of the momentum equation eliminates the pressure gradient term, yielding:
\begin{align}
\frac{\partial w}{\partial t} + \mathbf{u} \cdot \nabla w = \nu \nabla^2 w + \mathbf{f}_w,
\end{align}
where:
* $\mathbf{u} \cdot \nabla w$ represents the nonlinear advection of vorticity,
* $\nu \nabla^2 w$ accounts for viscous diffusion,
* $\mathbf{f}_w$ is the vorticity-specific forcing term.

The term $\mathbf{u} \cdot \nabla w$ can be expanded using the velocity components as:
\begin{align}
\mathbf{u} \cdot \nabla w = u_x \frac{\partial w}{\partial x} + u_y \frac{\partial w}{\partial y}.
\end{align}

By substituting $u_x$ and $u_y$ in terms of $\psi$, the nonlinear advection term is rewritten as the **Jacobian determinant**:
\begin{align}
J(\psi, w) = \frac{\partial \psi}{\partial x} \frac{\partial w}{\partial y} - \frac{\partial \psi}{\partial y} \frac{\partial w}{\partial x}.
\end{align}
Thus, the vorticity transport equation becomes:
\begin{align}
\frac{\partial w}{\partial t} - J(\psi, w) = \nu \nabla^2 w + \mathbf{f}_w.
\end{align}

+++

### Incorporating Additional Physical Effects

**Ekman damping** models frictional effects caused by the interaction of the fluid with a boundary layer.
It acts as a large-scale energy sink and is proportional to the vorticity:
\begin{align}
-\mu w,
\end{align}
where $\mu$ is the Ekman coefficient.
Including this term in the vorticity transport equation gives:
\begin{align}
\frac{\partial w}{\partial t} - J(\psi, w) = \nu \nabla^2 w - \mu w + \mathbf{f}_w.
\end{align}

Ekman damping is particularly relevant in geophysical systems, where it represents energy dissipation due to the Earth's surface or ocean floors.

The **$\beta$-plane approximation** models the variation of the Coriolis parameter $f$ with latitude.
In the vorticity equation, this introduces a term proportional to the northward velocity component $u_y$:
\begin{align}
\beta u_y,
\end{align}
where $\beta$ is the linear expansion coefficient of the Coriolis parameter.
Including this term in the vorticity transport equation gives:
\begin{align}
\frac{\partial w}{\partial t} - J(\psi, w) + \beta u_y = \nu \nabla^2 w - \mu w + \mathbf{f}_w.
\end{align}

The beta-plane term is crucial for studying large-scale atmospheric and oceanic dynamics, as it leads to phenomena such as Rossby waves and geostrophic turbulence.

### Advantages of the Vorticity-Streamfunction Formulation

1. **Elimination of Pressure:**
   The pressure term, which requires solving an additional Poisson equation in the velocity-pressure formulation, is completely removed in the vorticity-streamfunction approach.

2. **Reduced Number of Variables:**
   By working with $w$ and $\psi$, the system is reduced to a single scalar equation for $w$ coupled with the Poisson equation for $\psi$.

3. **Natural Compatibility with Spectral Methods:**
   The vorticity-streamfunction formulation lends itself well to spectral methods. Derivatives of $w$ and $\psi$ are straightforward to compute in spectral space, and the Poisson equation becomes an algebraic equation.

4. **Incorporation of Geophysical Effects:**
   Ekman damping and the beta-plane approximation extend the applicability of the formulation to real-world problems in atmospheric and oceanic sciences.

+++

## Spectral Representation of the Equations

The vorticity-streamfunction formulation simplifies the governing equations of 2D incompressible hydrodynamics, making them naturally suited to spectral methods.
By representing the vorticity and streamfunction in terms of Fourier series, we can efficiently compute derivatives and nonlinear terms in spectral space.
This section details the spectral representation of the vorticity transport equation and the Poisson equation, highlighting the mathematical transformations, strategies to handle aliasing errors, and implementation in Python.

+++

### Fourier Representation of Vorticity and Streamfunction

For a periodic domain of size $L_x \times L_y$, the vorticity $w(x, y, t)$ and streamfunction $\psi(x, y, t)$ can be expanded as Fourier series:
\begin{align}
w(x, y, t)    &= \sum_{k_x, k_y} \hat{w   }_{k_x, k_y}(t) e^{i (k_x x + k_y y)}, \\
\psi(x, y, t) &= \sum_{k_x, k_y} \hat{\psi}_{k_x, k_y}(t) e^{i (k_x x + k_y y)},
\end{align}
where:
* $\hat{w}_{k_x, k_y}$ and $\hat{\psi}_{k_x, k_y}$ are the Fourier coefficients of vorticity and streamfunction, respectively.
* $k_x$ and $k_y$ are the wavenumbers in the $x$- and $y$-directions, given by:
  \begin{align}
  k_x = \frac{2\pi n_x}{L_x}, \quad k_y = \frac{2\pi n_y}{L_y}, \quad n_x, n_y \in \mathbb{Z}.
  \end{align}

+++

#### Fourier Differentiation

In spectral space, derivatives with respect to $x$ and $y$ transform into simple multiplications by $ik_x$ and $ik_y$, respectively.
For example:
\begin{align}
\frac{\partial w}{\partial x} \rightarrow ik_x \hat{w}_{k_x, k_y}, \quad \frac{\partial^2 w}{\partial x^2} \rightarrow -k_x^2 \hat{w}_{k_x, k_y}.
\end{align}

This property makes spectral methods computationally efficient, as derivatives reduce to element-wise multiplications.

+++

### Poisson Equation in Spectral Space

The Poisson equation relates the vorticity $w$ and the streamfunction $\psi$:
\begin{align}
w = -\nabla^2 \psi.
\end{align}

In Fourier space, the Laplacian $\nabla^2$ becomes multiplication by $-k^2$, where $k^2 = k_x^2 + k_y^2$.
Thus, the Poisson equation transforms into:
\begin{align}
\hat{\psi}_{k_x, k_y} = \frac{\hat{w}_{k_x, k_y}}{k^2}, \quad k^2 \neq 0.
\end{align}

The $k^2 = 0$ mode corresponds to the mean value of $\psi$, which can be set to zero for flows with no net circulation.

+++

\begin{align}
u_y = -\frac{\partial \psi}{\partial x}.
\end{align}

+++

### Vorticity Transport Equation in Spectral Space

Recalling the vorticity transport equation in real space is:
\begin{align}
\frac{\partial w}{\partial t} - J(\psi, w) + \beta u_y = \nu \nabla^2 w - \mu w + \mathbf{f}_w.
\end{align}

In spectral space:
* The Laplacian term $\nabla^2 w$ transforms into $-k^2 \hat{w}_{k_x, k_y}$.
* The Ekman damping term $\mu w$ becomes $-\mu \hat{w}_{k_x, k_y}$.
* The Coriolis term $\beta u_y$ becomes $-\beta ik_x \hat{w}_{k_x, k_y}/k^2$.
* The nonlinear term $J(\psi, w)$ is evaluated in real space and transformed back to spectral space using the Fast Fourier Transform (FFT).

The equation in spectral space becomes:
\begin{align}
\frac{\partial}{\partial t}\hat{w}_{k_x, k_y} = \widehat{J(\psi, w)} + \beta \frac{ik_x \hat{w}_{k_x, k_y}}{k^2} - \nu k^2 \hat{w}_{k_x, k_y} - \mu \hat{w}_{k_x, k_y} + \hat{f}_{w k_x, k_y}.
\end{align}

+++

### Handling Aliasing Errors

Nonlinear terms, such as $J(\psi, w)$, involve products in real space that translate into convolutions in spectral space.
These convolutions can introduce spurious interactions between modes, known as **aliasing errors**, due to the finite resolution of the grid.

The **2/3 rule** is a widely used de-aliasing technique that truncates Fourier modes beyond $2/3$ of the maximum wavenumber.
For a grid with $N$ points in each dimension:
* Retain modes for $|k_x|, |k_y| \leq N/3$.
* Set all other modes to zero.

The 2/3 rule ensures that spurious contributions from nonlinear interactions fall outside the resolved spectral range.

+++

### Spectral-Galerkin vs. Pseudo-Spectral Methods

The **spectral-Galerkin method** projects the governing equations onto the basis functions of the spectral expansion.
This ensures the residual is orthogonal to the retained modes.

The **pseudo-spectral method** evaluates nonlinear terms in real space and transforms them back to spectral space using FFT.
While computationally efficient, it often requires de-aliasing or explicit viscosity to control aliasing errors.
This method balances speed and accuracy, making it popular for practical simulations.

```{code-cell} ipython3

```
