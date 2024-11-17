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

+++

### Stability Analysis Using the Linear Test Equation

To analyze the stability of numerical integrators, we commonly use a linear test equation from last lecture:
\begin{align}
\frac{dx}{dt} = \lambda x
\end{align}
where $\lambda \in \mathbb{C}$.
The exact solution is:
\begin{align}
x(t) = x_0 e^{\lambda t}.
\end{align}

+++

Consider a one-step numerical method applied to the test equation.
The update can generally be expressed as:
\begin{align}
x_{n+1} = R(z) x_n
\end{align}
where $R(z)$ is the amplification factor and $z = \lambda \Delta t$ is the stability parameter.

For the numerical method to be stable, the magnitude of the amplification factor must satisfy:
\begin{align}
|R(z)| \leq 1
\end{align}
This **stability condition** ensures that errors do not grow exponentially with each step.
