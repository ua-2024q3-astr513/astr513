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

Symbolic differentiation computes the derivative of a function expressed symbolically by applying calculus rules directly to its mathematical expression.
Unlike numerical methods that approximate derivatives at specific points, symbolic differentiation yields exact analytical expressions, making it valuable for theoretical analyses and precise computations.

## Algorithmic Approach

The general algorithm for symbolic differentiation involves:

1. Parsing the Expression: Represent the function as an expression tree, where each node corresponds to an operation (e.g., addition, multiplication) or operand (e.g., variable, constant).
2. Applying Differentiation Rules: Recursively apply differentiation rules to each node in the expression tree, including the chain rule.
3. Simplifying the Result: After applying the differentiation rules, simplify the resulting expression to make it more readable and computationally efficient.

Consider the function $f(x) = x^2 \sin(x) + e^{2x}$.
To compute $f'(x)$, a symbolic differentiation system would:
1. Differentiate $x^2 \sin(x)$ using the product rule:
   \begin{align}
   \frac{d}{dx}[x^2 \sin(x)] = x^2 \cos(x) + 2 x \sin(x)
   \end{align}
2. Differentiate $e^{2x}$ using the chain rule:
   \begin{align}
   \frac{d}{dx}[e^{2x}] = 2 e^{2x}
   \end{align}
3. Combine the results:
   \begin{align}
   f'(x) = x^2 \cos(x) + 2 x \sin(x) + 2 e^{2x}
   \end{align}

## Symbolic Computation with SymPy

`SymPy` is an open-source Python library for symbolic mathematics.
It allows for symbolic differentiation and manipulation of mathematical expressions.

Using `SymPy` to Compute $f'(x)$:

```{code-cell} ipython3
import sympy as sp

x = sp.symbols('x')
f = x**2 * sp.sin(x) + sp.exp(2 * x)
f_prime = sp.diff(f, x)
f_prime_simplified = sp.simplify(f_prime)

print("Derivative of f(x):")
sp.pprint(f_prime_simplified)
```

### Advantages of Symbolic Differentiation

Symbolic differentiation provides exact results by producing precise analytical expressions without approximation errors.
This exactness is crucial in theoretical work where precise solutions are necessary.
Additionally, symbolic derivatives are valid over continuous ranges of variables, offering general applicability in analyses.
They facilitate further analytical processes by simplifying tasks such as solving differential equations and optimizing functions, as the exact expressions can be manipulated algebraically to gain deeper insights into the behavior of the system under study.

+++

### Limitations

Symbolic differentiation has limitations that can impact its practicality.
One challenge is the potential for expression growth of the derivative expressions as the original function's complexity increases.
This growth can make the expressions difficult to interpret and computationally intensive to process.
The computational complexity associated with differentiating complex functions can require substantial resources and time, especially for high-dimensional systems or functions involving intricate compositions.
Furthermore, symbolic differentiation is not suitable for functions without explicit symbolic forms, such as those defined by experimental data, simulations, or complex numerical algorithms.
In such cases, alternative methods like numerical differentiation or automatic differentiation are more appropriate.

+++

### Software Tools

Several software tools facilitate symbolic differentiation by automating the application of calculus rules to mathematical expressions:

* [`SymPy`](https://www.sympy.org/):
  An open-source Python library that provides capabilities for symbolic differentiation, integration, and equation solving within the Python ecosystem.
* [`Mathematica`](https://www.wolfram.com/mathematica/):
  A computational software developed by Wolfram Research, offering extensive symbolic computation features used widely in academia and industry.
* [`Maple`](https://www.maplesoft.com/):
  A software package designed for symbolic and numeric computing, providing powerful tools for mathematical analysis.
* [`Maxima`](https://maxima.sourceforge.io/):
  An open-source computer algebra system specializing in symbolic manipulation, accessible for users seeking free alternatives.

+++

## Numerical Differentiation

Numerical differentiation estimates the derivative of a function using discrete data points, providing approximate values where analytical derivatives are difficult or impossible to obtain.
Unlike symbolic differentiation, which yields exact expressions, numerical methods offer flexibility in handling complex, empirical, or high-dimensional functions by leveraging computational algorithms to approximate derivatives at specific points.

+++

### Finite Difference Methods

Finite difference methods are fundamental techniques in numerical differentiation, estimating derivatives by evaluating the function at specific points and computing the ratio of differences.
These methods are essential when analytical derivatives are difficult or impossible to obtain, particularly for complex or empirical functions encountered in scientific and engineering applications.

The core idea behind finite difference methods is to approximate the derivative $f'(x)$ by evaluating the function $f(x)$ at selected points around $x$.
The three most basic finite difference approximations are forward difference, backward difference, and central difference.

The **forward difference approximation** estimates the first derivative at a point $x$ by using the function values at $x$ and $x + h$, where $h$ is a small step size:
\begin{align}
f'(x) \approx \frac{f(x + h) - f(x)}{h}.
\end{align}
This method is straightforward to implement and requires only one additional function evaluation beyond the point of interest.
However, its accuracy is limited by a truncation error of order $\mathcal{O}(h)$.
As $h$ decreases, the approximation becomes more accurate, but excessively small values of $h$ can lead to significant round-off errors due to the limitations of floating-point arithmetic.

Similarly, the **backward difference approximation** estimates the derivative using the function values at $x$ and $x - h$:
\begin{align}
f'(x) \approx \frac{f(x) - f(x - h)}{h}.
\end{align}
Like the forward difference, the backward difference method has the same truncation error of $\mathcal{O}(h)$.
It is particularly useful in situations where function evaluations at points greater than $x$ are not available or are computationally expensive.

The central difference approximation provides a more accurate estimate by averaging the forward and backward differences:
\begin{align}
f'(x) \approx \frac{1}{2}\left[\frac{f(x + h) - f(x)}{h} + \frac{f(x) - f(x - h)}{h}\right] = \frac{f(x + h) - f(x - h)}{2h}.
\end{align}
This method reduces the truncation error to $\mathcal{O}(h^2)$, making it  more accurate for smooth functions.
The central difference requires function evaluations at both $x + h$ and $x - h$, effectively doubling the number of required computations compared to the forward or backward methods.
Nevertheless, the enhanced accuracy often justifies the additional computational effort, especially in applications demanding high precision.

+++

### Error Analysis

Finite difference methods involve a trade-off between truncation error and round-off error.
The truncation error arises from approximating the derivative using discrete differences, while the round-off error is due to the finite precision of floating-point arithmetic used in computations.

For the forward and backward difference methods, the truncation error is proportional to $h$, meaning that decreasing $h$ improves the approximation's accuracy linearly.
In contrast, the central difference method’s truncation error decreases quadratically with $h$, offering better accuracy for smaller step sizes.

However, reducing $h$ too much can amplify round-off errors, as the difference $f(x + h) - f(x)$ becomes dominated by floating-point precision limitations.
Therefore, selecting an optimal step size $h$ is crucial.
Typically, $h$ is chosen to balance the minimization of both truncation and round-off errors, often on the order of the square root of machine epsilon (e.g., $h \approx \sqrt{\epsilon}$), where machine epsilon represents the smallest difference recognizable by the floating-point system.

+++

# Sample Implementation

To illustrate finite difference methods, consider the following Python implementation of the forward, backward, and central difference approximations:

```{code-cell} ipython3
def forward_difference(f, x, h):
    return (f(x + h) - f(x)) / h

def backward_difference(f, x, h):
    return (f(x) - f(x - h)) / h

def central_difference(f, x, h):
    return (f(x + h) - f(x - h)) / (2 * h)
```

```{code-cell} ipython3
import numpy as np

def f(x):
    return np.sin(x)

# Point of interest
def errs(x0 = np.pi / 4):

    # True derivative
    df0 = np.cos(x0)

    # Step sizes
    hs = np.logspace(0, -15, 31)

    errs_forward  = []
    errs_backward = []
    errs_central  = []

    for h in hs:
        df_forward  = forward_difference (f, x0, h)
        df_backward = backward_difference(f, x0, h)
        df_central  = central_difference (f, x0, h)
    
        errs_forward .append(abs(df_forward  - df0))
        errs_backward.append(abs(df_backward - df0))
        errs_central .append(abs(df_central  - df0))

    return hs, errs_forward, errs_backward, errs_central
```

```{code-cell} ipython3
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1,3, figsize=(12, 4), sharey=True)

for i, x0 in enumerate([0, np.pi/4, np.pi/2]):

    hs, errs_forward, errs_backward, errs_central = errs(x0)

    axes[i].loglog(hs, errs_forward,  'o-',  label='Forward Difference')
    axes[i].loglog(hs, errs_backward, 's--', label='Backward Difference')
    axes[i].loglog(hs, errs_central,  '^:',  label='Central Difference')
    axes[i].loglog(hs, hs,    lw=1)
    axes[i].loglog(hs, hs**2, lw=1)
    axes[i].set_xlim(1e1, 1e-16)
    axes[i].set_ylim(1e-17, 1e0)
    axes[i].set_xlabel('Step size h')
    axes[i].grid(True, which="both", ls="--")
    
axes[0].set_ylabel('Absolute Error')
axes[0].legend()
```

Why do the convergence rates do not behave as expected?

+++ {"jp-MarkdownHeadingCollapsed": true}

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
