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

# Root Finding and Optimization Techniques

Root finding and optimization are core numerical techniques that enable us to solve complex equations and optimize functions in fields where analytical solutions are often impossible.
Root finding aims to determine values for which a function $f(x) = 0$, and finds application across engineering, physics, and finance—whether calculating stresses in materials, energy levels in quantum mechanics, or rates of return in investments.
Optimization seeks to find the minimum or maximum of a function and is especially crucial in machine learning, where minimizing loss functions directly affects model performance.
The two concepts intersect in gradient-based optimization, where finding the roots of a gradient helps locate stationary points and optimize complex models.

+++

## Root Finding Techniques

### Bisection Method

The Bisection Method is a simple and robust root-finding algorithm that relies on the Intermediate Value Theorem.
The theorem states that if $f(x)$ is a continuous function on an interval $[a, b]$ and $f(a)$ and $f(b)$ have opposite signs, then there exists at least one root in the interval $(a, b)$ where $f(x) = 0$.
We already implemented a similar algorithm in a [previous lecture](interpolate.md).

```{code-cell} ipython3
def bisection_search(xs, target):
    l, h = 0, len(xs) - 1
    while h - l > 1:
        m = (h + l) // 2
        if target >= xs[m]:
            l = m
        else:
            h = m
    return l # returns index of the closest value less than or equal to target
```

The main difference is that we no longer have a finite set of sampling points.

```{code-cell} ipython3
def bisection(f, l, h, tol=1e-6):
    if f(l) * f(h) >= 0:
        raise ValueError("f(a) and f(b) must have opposite signs")
    while h - l > 2*tol:
        m = (l + h) / 2
        if f(m) == 0:
            return m  # c is the root
        elif f(l) * f(m) > 0:
            l = m
        else:
            h = m
    return (l + m) / 2
```

Example usage:

```{code-cell} ipython3
def f(x):
    return x**3 - x - 2

root = bisection(f, 1, 2)
print("Approximate root:")
print("  x0  = ",   root )
print("f(x0) = ", f(root))
```

### Newton-Raphson Method

### Secant Method

### Van Wijngaarden-Dekker-Brent Method

### Multidimensional Root Finding

## Optimization Methods

### Gradient Descent Methods

### Stochastic Gradient Descent (SGD)

### Momentum and Adaptive Methods

### Conjugate Gradient Methods

### Lagrange Multipliers and KKT Conditions

### Simulated Annealing

## Optimization in Machine Learning

### Role of Optimization in Training Models

### Challenges in Deep Learning

### Hyperparameter Optimization

## Connecting Root Finding and Optimization

## Conclusion

```{code-cell} ipython3

```
