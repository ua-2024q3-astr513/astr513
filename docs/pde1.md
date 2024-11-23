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

### From Boltzmann to Navier-Stokes

While the integral-to-differential derivation provides a macroscopic view of fluid dynamics, an alternative microscopic approach begins with the Boltzmann Equation and employs the Moment Method to derive the Navier-Stokes equations.
This method bridges the gap between kinetic theory and continuum mechanics, offering deeper insights into the assumptions underlying fluid behavior.

+++

#### Boltzmann Equation Overview

The Boltzmann Equation describes the statistical behavior of a dilute gas out of equilibrium.
It is given by:
\begin{align}
\frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla f + \frac{\mathbf{F}}{m} \cdot \frac{\partial f}{\partial \mathbf{v}} = \left( \frac{\partial f}{\partial t} \right)_{\text{coll}}
\end{align}
where:
* $f = f(\mathbf{x}, \mathbf{v}, t)$ is the distribution function,
* $\mathbf{F}$ is an external force,
* $m$ is the particle mass,
* The right-hand side represents the collision term.

+++

#### Moment Method

To derive macroscopic equations from the Boltzmann Equation, we take moments by multiplying the equation by powers of velocity $\mathbf{v}$ and integrating over all velocities.

1. **Zeroth Moment (Continuity Equation):**
   \begin{align}
   \int \left( \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla f \right) d\mathbf{v} = 0 \quad \Rightarrow \quad \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0
   \end{align}

2. **First Moment (Momentum Equation):**
   \begin{align}
   \int \mathbf{v} \left( \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla f \right) d\mathbf{v} = \int \mathbf{F} f \, d\mathbf{v} - \int \mathbf{v} \left( \frac{\partial f}{\partial t} \right)_{\text{coll}} d\mathbf{v}
   \end{align}
   Simplifying leads to:
   \begin{align}
   \frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u}) = -\nabla p + \nabla \cdot \mathbf{\Pi} + \rho \mathbf{f}
   \end{align}

3. **Second Moment (Energy Equation):**
   \begin{align}
   \int \frac{1}{2} m v^2 \left( \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla f \right) d\mathbf{v} = \int \frac{1}{2} m v^2 \mathbf{F} \cdot \nabla_{\mathbf{v}} f \, d\mathbf{v} - \int \frac{1}{2} m v^2 \left( \frac{\partial f}{\partial t} \right)_{\text{coll}} d\mathbf{v}
   \end{align}
   This results in:
   \begin{align}
   \frac{\partial (\rho e)}{\partial t} + \nabla \cdot (\rho e \mathbf{u}) = \nabla \cdot \mathbf{q} + \Phi + \rho \mathbf{f} \cdot \mathbf{u}
   \end{align}

These moments collectively yield the **Navier-Stokes Equations**, providing a comprehensive description of fluid motion from a kinetic theory perspective.

+++

### Assumptions & Physical Meaning

Several key assumptions underpin the derivation of the Navier-Stokes equations, simplifying the complex interactions within a fluid to make the equations tractable:

1. **Continuum Hypothesis:**
   
   Assumes that fluid properties such as density and velocity are continuously distributed in space.
   This is valid when the mean free path of particles is much smaller than the characteristic length scales of interest.

2. **Newtonian Fluid:**
   
   Assumes that the stress tensor $\mathbf{\Pi}$ is linearly related to the strain rate tensor.
   This implies that the fluid's viscosity is constant and that the relationship between shear stress and strain rate is proportional.

4. **Local Thermodynamic Equilibrium:**
   
   Assumes that the distribution function $f$ is close to the Maxwell-Boltzmann distribution, allowing for the expansion of $f$ in terms of small deviations.
   This ensures that macroscopic quantities like temperature and pressure are well-defined locally.

5. **Negligible External Forces:**
   
   Often, body forces $\mathbf{f}$ such as gravity are considered constant or negligible compared to other forces.
   In cases where external forces are significant, they are incorporated explicitly into the momentum and energy equations.

6. **Incompressible Flow:**
   
   For many practical applications, especially at low Mach numbers, the fluid is assumed to be incompressible ($\nabla \cdot \mathbf{u} = 0$).
   This simplifies the continuity equation and decouples the pressure from the density.

+++

#### Viscosity, Pressure Gradients, and Other Key Terms

- **Pressure Gradient ($-\nabla p$):**
  
  Represents the force per unit volume exerted by pressure differences within the fluid.
  It drives the fluid from regions of high pressure to low pressure, influencing the acceleration and direction of flow.

- **Viscous Terms ($\mu \nabla^2 \mathbf{u}$):**
  
  Account for internal friction due to the fluid's viscosity, which resists deformation and dissipates kinetic energy into thermal energy.
  Viscosity is a measure of a fluid's resistance to shear or flow.

- **Body Forces ($\rho \mathbf{f}$):**
  
  Represent external forces acting on the fluid, such as gravity or electromagnetic forces.
  These forces can influence the overall motion and stability of the fluid flow.

These terms collectively describe the balance between inertial forces, pressure forces, viscous forces, and external influences, shaping the fluid's behavior under various conditions.

+++

#### Limiting Cases and Their Physical Interpretations

Analyzing the Navier-Stokes equations under different assumptions and limiting cases provides deeper insights into fluid behavior:

1. **Inviscid Flow (Euler Equations):**
   
   By setting viscosity $\mu = 0$, the Navier-Stokes equations reduce to the Euler equations:   
   \begin{align}
   \rho \left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} \right) = -\nabla p + \rho \mathbf{f}
   \end{align}
   These equations describe ideal, non-viscous fluid flow and are applicable in scenarios where viscous effects are negligible, such as high Reynolds number flows.

2. **Steady-State Flow:**
   
   Assuming no temporal changes ($\frac{\partial}{\partial t} = 0$), the equations describe steady-state conditions where fluid properties remain constant over time.
   This simplification is useful in analyzing flows that have reached equilibrium.

4. **Laminar vs. Turbulent Flow:**
   
   * **Laminar Flow:** Occurs at low Reynolds numbers where viscous forces dominate, resulting in smooth, orderly fluid motion.
   * **Turbulent Flow:** Arises at high Reynolds numbers where inertial forces prevail, leading to chaotic and eddying fluid motion.

Understanding these limiting cases aids in selecting appropriate models and computational approaches for different fluid flow regimes.

+++

### Significance in Fluid Dynamics

The derived fluid dynamics equations, particularly the Navier-Stokes equations, are foundational in describing and predicting fluid behavior across a myriad of applications.
Their significance extends from everyday engineering problems to the most complex astrophysical phenomena.

+++

#### Role of PDEs in Governing Real-World Flows

Partial Differential Equations (PDEs) like the Navier-Stokes equations encapsulate the fundamental principles governing fluid motion—conservation of mass, momentum, and energy.
These equations enable the simulation and analysis of fluid flows in various contexts:

* **Atmospheric Sciences:**
  
  PDEs model weather patterns, climate dynamics, and the movement of air masses, aiding in weather forecasting and climate modeling.

* **Oceanography:**
  
  They describe ocean currents, wave dynamics, and the mixing of different water masses, contributing to the understanding of marine ecosystems and global climate.

* **Engineering Applications:**
  
  In aerospace engineering, PDEs are used to design and optimize aircraft and spacecraft by simulating airflow around structures.
  In civil engineering, they aid in modeling the flow of water in pipes and channels.

* **Astrophysics:**
  
  PDEs are indispensable for modeling the behavior of astrophysical fluids under extreme conditions, such as those found in stellar interiors, accretion disks, and interstellar gas clouds.

+++

#### Application in Astrophysical Fluid Systems

Astrophysical systems often involve complex fluid dynamics under extreme conditions of temperature, pressure, and gravity.
PDEs are essential in modeling these systems, providing insights into their structure, evolution, and behavior.

1. **Stellar Interiors:**
   
   The interiors of stars involve convective and radiative energy transport, described by the Navier-Stokes equations coupled with the heat equation.
   These models help in understanding stellar structure, energy generation, and evolution.

2. **Accretion Disks:**
   
   Around compact objects like black holes and neutron stars, accretion disks form as matter spirals inward.
   PDEs model the fluid dynamics, angular momentum transport, and radiation processes within these disks, crucial for explaining phenomena like quasars and X-ray binaries.

3. **Interstellar Medium (ISM):**
   
   The ISM consists of gas and dust influenced by gravitational, magnetic, and thermal forces.
   Magnetohydrodynamic (MHD) equations, a set of PDEs, describe the behavior of conducting fluids in the ISM, aiding in the study of star formation, supernova remnants, and galactic dynamics.

4. **Gravitational Waves:**
   
   The propagation of gravitational waves, ripples in spacetime caused by massive astrophysical events, is described by the Einstein Field Equations, which are highly complex nonlinear PDEs.
   Numerical solutions to these equations are essential for predicting gravitational wave signatures detected by observatories like LIGO and Virgo.

5. **Supernova Explosions:**
   
   The explosive death of massive stars involves shock waves, turbulence, and nuclear reactions, all governed by fluid dynamics PDEs.
   Modeling these processes is crucial for understanding the distribution of heavy elements and the dynamics of supernova remnants.

+++

## Classification of Partial Differential Equations (PDEs)

Partial Differential Equations (PDEs) are integral to modeling a vast array of physical phenomena, from the diffusion of heat to the propagation of waves.
Understanding the classification of PDEs is essential for selecting appropriate solution methods and interpreting the behavior of the systems they describe.
PDEs are primarily classified into three categories: **Elliptic**, **Parabolic**, and **Hyperbolic**.
Each class possesses distinct mathematical properties and governs different types of physical phenomena.
This section delves into the nature, applications, and astrophysical relevance of each PDE type.

+++

### Elliptic PDEs

**Definition and Characteristics:**

Elliptic PDEs are characterized by the absence of real characteristics, implying that the solution at any given point depends on the values of the solution in an entire surrounding region.
These equations typically describe steady-state situations where there is no time dependence.

The general second-order linear elliptic PDE in two variables can be written as:
\begin{align}
A \frac{\partial^2 u}{\partial x^2} + 2B \frac{\partial^2 u}{\partial x \partial y} + C \frac{\partial^2 u}{\partial y^2} + \text{lower order terms} = 0
\end{align}
where the discriminant $B^2 - AC < 0$, ensuring ellipticity.

**Example: Laplace's Equation:**

One of the most fundamental elliptic PDEs is Laplace's equation, given by:
\begin{align}
\nabla^2 u = 0
\end{align}
where $\nabla^2$ is the Laplacian operator.
Laplace's equation arises in various contexts, including electrostatics, gravitational potentials, and incompressible fluid flow.

**Physical Intuition and Applications:**

Elliptic PDEs model phenomena where the state of the system is in equilibrium.
Solutions to elliptic equations are generally smooth and exhibit no abrupt changes, reflecting the steady-state nature of the systems they describe.

**Astrophysical Applications:**

In astrophysics, Laplace's equation is used to model gravitational potentials in regions devoid of mass.
For instance, the gravitational potential outside a spherically symmetric mass distribution can be described by Laplace's equation, facilitating the study of stellar and planetary structures.

+++

#### Parabolic PDEs

**Definition and Characteristics:**

Parabolic PDEs are associated with diffusion-like processes and typically involve one temporal dimension and multiple spatial dimensions.
They describe how a system evolves over time towards equilibrium.

The general form of a second-order linear parabolic PDE in two variables is:
\begin{align}
\frac{\partial u}{\partial t} = A \frac{\partial^2 u}{\partial x^2} + B \frac{\partial u}{\partial x} + C u + \text{source terms}
\end{align}
with the discriminant $B^2 - 4AC = 0$.

**Example: Heat Equation:**

The quintessential parabolic PDE is the heat equation:
\begin{align}
\frac{\partial u}{\partial t} = \alpha \nabla^2 u
\end{align}
where $u$ represents temperature, and $\alpha$ is the thermal diffusivity of the material.

**Physical Intuition and Applications:**

Parabolic PDEs model processes where diffusion plays a central role, leading to the smoothing out of gradients over time.
Solutions to parabolic equations typically show how an initial distribution evolves, gradually approaching a steady state.

**Astrophysical Applications:**

In astrophysics, parabolic PDEs are employed to model thermal conduction within stars.
For example, the transport of heat from the core to the outer layers of a star can be described by the heat equation, aiding in the understanding of stellar evolution and energy distribution.

+++

#### Hyperbolic PDEs

**Definition and Characteristics:**

Hyperbolic PDEs are associated with wave propagation and are characterized by the existence of real characteristics along which information travels.
These equations typically involve both temporal and spatial derivatives and can model dynamic, time-dependent phenomena.

The general form of a second-order linear hyperbolic PDE in two variables is:
\begin{align}
A \frac{\partial^2 u}{\partial t^2} + B \frac{\partial^2 u}{\partial t \partial x} + C \frac{\partial^2 u}{\partial x^2} + \text{lower order terms} = 0
\end{align}
where the discriminant $B^2 - 4AC > 0$, ensuring hyperbolicity.

**Example – Wave Equation:**

A fundamental example of a hyperbolic PDE is the wave equation:
\begin{align}
\frac{\partial^2 u}{\partial t^2} = c^2 \nabla^2 u
\end{align}
where $u$ represents the wave function, and $c$ is the wave speed.

**Physical Intuition and Applications:**

Hyperbolic PDEs model systems where disturbances propagate as waves, maintaining their shape over time.
Solutions to hyperbolic equations exhibit finite propagation speeds and can describe phenomena such as sound waves, electromagnetic waves, and shock waves.

**Astrophysical Applications:**

In the astrophysical context, hyperbolic PDEs are essential for modeling wave phenomena in stellar oscillations and the propagation of shock waves in supernova explosions.
Additionally, the Einstein Field Equations, which describe gravitational waves in general relativity, are a set of highly nonlinear hyperbolic PDEs.

+++

### Physical Intuition Behind Each PDE Type

Understanding the physical intuition behind each class of PDEs aids in selecting appropriate models and solution techniques for different phenomena:

* **Elliptic PDEs** describe equilibrium states where the system is balanced, and there is no inherent time evolution. They are often used to determine steady-state distributions of physical quantities.

* **Parabolic PDEs** capture the time-dependent evolution of systems towards equilibrium through diffusive processes. They are ideal for modeling heat conduction, mass diffusion, and similar phenomena.

* **Hyperbolic PDEs** represent dynamic systems where disturbances propagate as waves. They are suited for modeling scenarios involving wave motion, vibrations, and other propagating phenomena.

+++

### Examples of Physical Phenomena

* **Elliptic PDEs:**
  * **Electrostatic Potential:** Determined by Laplace's or Poisson's equation in regions without or with charge distributions, respectively.
  * **Steady-State Fluid Flow:** Described by the incompressible Navier-Stokes equations under steady conditions.
  
* **Parabolic PDEs:**
  * **Heat Conduction:** Governed by the heat equation, modeling temperature distribution over time.
  * **Diffusion Processes:** Including the spread of pollutants in the atmosphere or solutes in a solvent.
  
* **Hyperbolic PDEs:**
  * **Sound Waves:** Described by the acoustic wave equation, modeling pressure variations in a medium.
  * **Electromagnetic Waves:** Governed by Maxwell's equations, which can be expressed as hyperbolic PDEs.
  * **Gravitational Waves:** Modeled by the Einstein Field Equations in the context of general relativity.

+++

## Non-Dimensionalization and Key Dimensionless Numbers

Non-dimensionalization is a fundamental technique in the analysis of partial differential equations (PDEs) and fluid dynamics.
By converting variables and equations into dimensionless forms, researchers and engineers can simplify complex systems, identify dominant physical effects, and generalize solutions across different scales.
This section explores the purpose of non-dimensionalization, introduces key dimensionless numbers—Reynolds Number (Re), Mach Number (Ma), and Prandtl Number (Pr)—and discusses their applications and significance in various physical systems.

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

+++ {"jp-MarkdownHeadingCollapsed": true}

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
