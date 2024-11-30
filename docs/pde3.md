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
