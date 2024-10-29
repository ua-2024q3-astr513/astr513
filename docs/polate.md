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

## Preliminaries: Searching an Ordered Table

## Polynomial Interpolation and Extrapolation

## Cubic Spline Interpolation

## Rational Function Interpolation and Extrapolation

## Coefficients of the Interpolating Polynomial

## Interpolation on a Grid in Multidimensions

## Interpolation on Scattered Data in Multidimensions

## Laplace Interpolation

## Conclusion and Discussion
