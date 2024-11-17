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

# Integration of Ordinary Differential Equations. II

+++

## Numerical Stability of Integrators

Numerical Stability in the context of ODE solvers refers to the ability of a numerical method to control the growth of errors introduced during the iterative process of approximation.
A method is stable if the numerical errors do not amplify uncontrollably as computations proceed.
* Definition: A numerical method is stable for a given problem if the errors (from truncation or round-off) do not grow exponentially with time steps.
* Importance: Stability ensures that the numerical solution behaves similarly to the true solution, especially over long integration intervals.

Stability is a different concept than accuracy.
* Stability: Concerns the boundedness of errors over time.
* Accuracy: Pertains to how closely the numerical solution approximates the exact solution.

A method can be stable but not necessarily accurate, and vice versa.
However, both properties are requried for reliable numerical solutions.
