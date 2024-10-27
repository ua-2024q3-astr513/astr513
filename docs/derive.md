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

# Numerical and Automatic Derivatives

Derivatives are fundamental in mathematical modeling, providing essential insights into the behavior of physical systems by quantifying rates of change.
In fields such as computational physics, engineering, and machine learning, the efficient and accurate computation of derivatives is crucial for simulations, optimizations, and analyses.
Traditional analytical methods for finding derivatives may become intractable for complex or nonlinear functions commonly encountered in real-world applications
Consequently, alternative techniques have emerged as indispensable tools in scientific computing.

The derivative of a real-valued function $f(x)$ at a point $x = a$ is defined as the limit:
\begin{align}
f'(a) = \lim_{h \to 0} \frac{f(a + h) - f(a)}{h}.
\end{align}
This limit, if it exists, represents the slope of the tangent line to the curve $y = f(x)$ at $x = a$. The derivative function $f'(x)$ provides the rate of change of $f$ at any point within its domain where the derivative exists.

Several fundamental rules facilitating the computation of derivatives are taught in undergraduate calculus courses.
Among them, the most important one is the chain rule.
It states that, for $f(x) = g(h(x))$, its derivative is given by
\begin{align}
f'(x) = g'(h(x)) h'(x).
\end{align}
We will show that the chain rule is extremely important in modern numerical and automatic derivatives.

Methods for computing derivatives include symbolic differentiation, numerical approximation, and automatic differentiation.
Symbolic differentiation applies analytical rules directly to mathematical expressions, yielding exact derivative formulas.
Numerical methods, such as finite difference techniques, approximate derivatives using discrete data points and are straightforward to implement but may suffer from truncation and round-off errors.
Automatic differentiation bridges the gap by systematically applying the chain rule to compute exact derivatives up to machine precision without symbolic manipulation, making it efficient for complex functions and large-scale systems.

Understanding the principles, advantages, and limitations of these approaches allows for the selection of the most appropriate method for a given problem.
This lecture will introduce these techniques, providing a comprehensive overview of their theoretical foundations and practical implementations in computational contexts.

+++

## Symbolic Differentiation

+++

## Numerical Differentiation

### Finite Difference Methods

* Low-Order Finite Difference Formulas
  * Forward difference formula.
  * Backward difference formula.
  * Central difference formula.
  * Error analysis and truncation errors.

* High-Order Finite Difference Methods
  * Derivation of higher-order formulas.
  * Error reduction and convergence rates.

### Spectral Methods and Fourier Transform

* Connection between finite differences and Fourier Transform.
* Introduction to spectral differentiation.
* Advantages for smooth periodic functions.

### Complex Step Differentiation

* Introduction to Complex Step Method
* Concept and mathematical foundation.
* Elimination of subtractive cancellation errors.

+++

## Automatic Differentiation

### Introduction to Automatic Differentiation

* Motivation and need for AD.
* Differences between AD, symbolic differentiation, and numerical approximation.

### Dual Numbers and Forward Mode AD

* Understanding Dual Numbers
  * Definition and algebra of dual numbers.
  * How dual numbers enable forward mode AD.
* Connection to Complex Step Method
  * Mathematical similarities and practical differences.

### Reverse Mode AD and Backpropagation

* Concept of Reverse Mode AD
  * Computational graphs and the chain rule.
* Backpropagation in Neural Networks
  * Application in machine learning.

### Higher-Order Derivatives: Jacobians and Hessians

* Importance of higher-order derivatives
  * Applications in optimization algorithms.
  * Computing Jacobians and Hessians with AD
* Techniques and computational considerations.

### Limitations and Challenges of AD

* Handling non-differentiable functions and discontinuities.
* Computational overhead and memory usage.
* Control flow and dynamic graphs.
* Best practices for efficient AD implementation.

+++

## Comparison of Differentiation Methods

### Numerical vs. Automatic vs. Symbolic Differentiation

* Strengths and weaknesses of each method.

### Choosing the right method for different applications.

+++

## Summary

### Open questions and discussion

### Suggestions for further reading and exploration.

```{code-cell} ipython3

```
