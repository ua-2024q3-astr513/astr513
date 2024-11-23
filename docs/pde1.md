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

## Introduction to Partial Differential Equations (PDEs) (10 minutes)

* What are PDEs?
  
* Definition and significance in modeling continuous systems.
  
* Contrast with Ordinary Differential Equations (ODEs).

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
  * Elliptic: Example – Laplace’s equation.
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
