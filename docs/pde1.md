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

## Introduction to PDEs (15 minutes)
- **What are PDEs?**
  - Definition and importance in modeling continuous systems.
  - Comparison with ordinary differential equations (ODEs).
  - Examples of real-world applications:
    - Fluid flow (Navier-Stokes).
    - Heat transfer (diffusion equation).
    - Wave propagation (wave equation).
- **Classification of PDEs**:
  - Elliptic: Laplace’s equation.
  - Parabolic: Heat equation.
  - Hyperbolic: Wave equation.
  - Discuss the physical intuition behind each type.

+++

## Derivation of Fluid Dynamics Equations (15 minutes)
- **From Boltzmann to Navier-Stokes**:
  - Briefly introduce the Boltzmann equation.
  - Use the moment method to derive:
    - Continuity equation.
    - Momentum equation.
    - Energy equation.
  - Emphasize the assumptions and physical meaning of terms (e.g., viscosity, pressure gradient).
- **Significance of PDEs in Fluid Dynamics**:
  - How these equations govern flows in real-world systems (e.g., airflow, ocean currents).

+++

## Non-Dimensionalization and Key Dimensionless Numbers (15 minutes)
- **Why Non-Dimensionalization?**
  - Simplifying equations for analysis.
  - Identifying dominant physical effects in specific regimes.
- **Key Dimensionless Numbers**:
  - Reynolds Number $\text{Re}$: Ratio of inertial to viscous forces.
  - Mach Number     $\text{Ma}$: Compressibility effects.
  - Prandtl Number  $\text{Pr}$: Momentum vs. thermal diffusivity.
  - Discuss their roles in determining the behavior of physical systems.

+++

## Numerical Techniques for Solving PDEs (30 minutes)
- **Finite Difference Methods**:
  - Forward Time Centered Space (FTCS):
    - Formulation and its instability.
  - von Neumann Stability Analysis:
    - Define the amplification factor and analyze stability.
- **Stabilized Methods**:
  - Lax Method: Adding numerical dissipation.
  - Courant-Friedrichs-Lewy (CFL) Condition: Ensuring stability.
- **Comparison of Methods**:
  - Highlight practical considerations for choosing a method (stability, accuracy, computational cost).

+++

## Practical Example: Heat Equation (15 minutes)
- **1D Heat Equation**:
  - Derive the equation $u_t = \alpha u_{xx}$ from conservation principles.
  - Solve analytically using separation of variables.
  - Discuss boundary and initial conditions.
- **Numerical Solution**:
  - Outline a basic finite difference implementation.
  - Discuss the impact of the CFL condition.

```{code-cell} ipython3

```
