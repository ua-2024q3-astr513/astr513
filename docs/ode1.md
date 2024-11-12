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

# Integration of Ordinary Differential Equations. I

In this lecture, we will explore numerical methods for solving ordinary differential equations (ODEs), focusing on initial value problems (IVPs).
We will begin by discussing the importance of integrating differential equations in physics and the limitations of analytical solutions.
We will then introduce Euler's method, both forward (explicit) and backward (implicit) forms.
We will examine how truncation errors and round-off errors affect these methods.
Additionally, we will demonstrate that the backward Euler method is more stable than the forward Euler method by analyzing their stability regions.

+++

## Introduction to Numerical ODEs

Ordinary differential equations (ODEs) are foundational to understanding and predicting the behavior of physical systems.
Many of the dynamic processes that occur in nature, from the motion of planets under gravitational forces to the oscillations of atoms in a crystal lattice, are described by differential equations.
For instance, Newton's second law, which states that the force on an object is equal to its mass times its acceleration, can be written as:
\begin{align}
F = m a = m \frac{d^2x}{dt^2}
\end{align}
This equation, along with others describing physical laws, typically involves derivatives of a function, representing rates of change, which relate the state of a system at one moment to its state at another.
In cases where physical fields, like the electromagnetic or gravitational fields, are described, partial differential equations (PDEs) are used since they involve rates of change with respect to multiple variables.
However, in this lecture, we focus on ODEs where functions depend on a single variable, typically time.

It is not always possible to solve ODEs analytically
Analytic solutions are often limited to simple or highly idealized cases.
Real-world systems tend to involve complex, nonlinear relationships that cannot be resolved through straightforward integration or algebraic manipulation.
For these cases, numerical methods provide a powerful alternative.

+++

### Problem Definition and Types of Differential Equations

Consider the following two types of differential equations, which differ in their dependence on variables:

1.  The first type has the form:
    \begin{align}
    \frac{dx}{dt} = f(t)
    \end{align}
    Here, the right-hand side (RHS) of the equation, $f(t)$, depends only on the independent variable $t$ and not on $x$, the dependent variable we are solving for.
2.	The second type is given by:
    \begin{align}
    \frac{dx}{dt} = f(x, t)
    \end{align}
    In this equation, the RHS depends on both $t$ and $x$, making it more complicated since $x$, the quantity we are trying to determine, appears in the expression that defines its rate of change.

+++

The first type of differential equation, where $f(t)$ depends only on $t$, can be solved by direct integration.
This is because we can integrate $f(t)$ with respect to $t$ to find $x(t)$:
\begin{align}
x(t) = x(t_0) + \int_{t_0}^t f(t') \, dt'
\end{align}
This approach is feasible when $f(t)$ is a function that can be integrated analytically.
In cases where this is not possible, numerical integration techniques, such as the trapezoidal rule or Simpson's rule, can be used to approximate the solution.
We covered these numerical integration techniques in an earlier lecture on "[Numerical Integration of Functions](integration.md)".

+++

The second type of differential equation, where $f(x, t)$ depends on both $x$ and $t$, is more complicated.
We cannot solve this type by direct integration because $x$, the function we are trying to find, is also part of the expression on the RHS. For example, attempting to directly integrate is not feasible since $x$ is unknown within the integral itself.
In other words, we do not know $x$ at intermediate points between $t_0$ and $t$, so we cannot compute the integral without first determining $x(t)$ at these points.

This is where numerical methods for ODEs become essential.
Instead of solving for $x(t)$ in one go, we use numerical methods to approximate the solution by advancing in small increments over the interval of interest.
This way, we iteratively approximate $x(t)$ at discrete points in time, allowing us to handle equations where direct integration is not possible.

Since many real-world systems fall into the second category of differential equations, where the dependence on  $x$ and $t$ is nonlinear, numerical methods are widely used in scientific computing.
These methods break down the problem into manageable steps and enable us to approximate solutions even when analytic solutions do not exist or are difficult to obtain.
Numerical methods not only provide practical solutions to ODEs but also allow for flexibility in modeling a wide range of phenomena, such as chaotic systems, non-linear oscillations, and biological systems with interacting populations.

```{code-cell} ipython3

```
