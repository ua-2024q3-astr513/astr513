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

## Introduction

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
