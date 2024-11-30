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
