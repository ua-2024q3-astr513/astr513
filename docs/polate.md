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

# Interpolation and Extrapolation

This lecture follows closely (2nd Edition in C and 3rd Edition in C++), Chapter 3 "Interpolation and Extrapolation".

+++

## Introduction

In scientific computing and machine learning, interpolation and extrapolation are essential tools for estimating function values at new data points based on known information.
In machine learning, all standard supervised learning tasks can be viewed as interpolation problems in high-dimensional space, where models predict outputs within the range of training data.
When attempting to predict outside this range, however, we enter the realm of extrapolation, often referred to as out-of-domain generalization in machine learning.
Extrapolation is challenging because models typically lack information beyond their training data, making reliable predictions difficult.

Interpolation methods include polynomial and rational function interpolation, as well as spline approaches.
Polynomial interpolation is versatile but prone to significant oscillations, especially at the edges of data (Runge’s phenomenon).
Rational functions, which use ratios of polynomials, can offer more stable estimates and handle asymptotic behavior better.
Spline interpolation, particularly cubic splines, is valued for its smoothness and continuity up to the second derivative, making it effective for applications requiring a smooth fit.

Extrapolation remains difficult, yet physics-informed machine learning (PIML) presents a promising avenue.
By embedding known physical laws, such as ordinary differential equations (ODEs), into models, PIML enables extrapolation that aligns with fundamental constraints, making it possible to extend predictions meaningfully beyond the observed data range.

Interpolation and function approximation are related but distinct tasks.
While interpolation estimates values at specified points within a given dataset, function approximation creates a simplified function to replace a more complex one.
In approximation, we can sample additional points as needed, whereas interpolation relies strictly on values at specific, fixed sampling points.
(See [Numerical Recipes](https://numerical.recipes/) Chapter 5 for function approximation.)

Interpolation also has limitations.
Pathological functions can defy even the most sophisticated interpolation schemes.
For example, consider a function that behaves smoothly except for a slight singularity at a certain point:
\begin{align}
f(x) = 3x^2 + \frac{1}{\pi^4}\ln\left[(\pi - x)^2\right] + 1
\end{align}
Interpolation based on values close to but not precisely at that singularity will likely produce an inaccurate result.

```{code-cell} ipython3
import numpy as np

def f(x):
    return 3 * x**2 + np.log((np.pi - x)**2) / np.pi**4 + 1

x1 = np.array([3.13, 3.14, 3.15, 3.16])
x2 = np.linspace(3.13, 3.16, 31)
x3 = np.linspace(3.13, 3.16, 301)
x4 = np.linspace(3.13, 3.16, 3001)
```

```{code-cell} ipython3
from matplotlib import pyplot as plt

plt.plot(x4, f(x4))
plt.plot(x3, f(x3), '--')
plt.plot(x2, f(x2), 'o:')
plt.plot(x1, f(x1), 'o-')
```

These cases highlight the importance of incorporating error estimates in interpolation routines.
Although no error estimate is foolproof, an effective interpolation method should still offer a reasonable assessment of its own accuracy within the presumption of smoothness.

+++

## Preliminaries: Searching an Ordered Table

In many interpolation tasks, especially with irregularly sampled data, the process begins with a critical first step: identifying the nearest points surrounding the target interpolation value.

Unlike regularly spaced data on a uniform grid, where adjacent points are easy to locate by simple indexing, randomly sampled or unevenly spaced data requires additional steps to find nearby values.
This searching step can be as computationally intensive as the interpolation itself, so efficient search methods are essential to maintain overall performance.

In Numerical Recipes, two primary methods are presented for this purpose: bisection and hunting.
Each is suited to different scenarios, depending on whether interpolation points tend to be close to one another or scattered randomly.

+++

## Polynomial Interpolation and Extrapolation

## Cubic Spline Interpolation

## Rational Function Interpolation and Extrapolation

## Coefficients of the Interpolating Polynomial

## Interpolation on a Grid in Multidimensions

## Interpolation on Scattered Data in Multidimensions

## Laplace Interpolation

## Conclusion and Discussion
