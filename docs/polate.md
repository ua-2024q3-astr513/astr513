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

### Linear Search

As a reference, we will implement a linear search:

```{code-cell} ipython3
def linear(xs, target):
    for l in range(len(xs)): # purposely use for-loop to avoid C optimization in numpy
        if xs[l] >= target:
            return l-1
```

```{code-cell} ipython3
import numpy as np

for _ in range(10):
    xs = np.sort(np.random.uniform(0, 100, 10))
    v  = np.random.uniform(min(xs), max(xs))
    i  = linear(xs, v)
    print(f'{xs[i]} <= {v} < {xs[i+1]}')
```

### Bisection Search

Bisection search is a reliable method that works by dividing the search interval in half with each step until the target value’s position is found.
Given a sorted array of $N$ data points, this method requires approximately $\log_2(N)$ steps to locate the closest point, making it efficient even for large datasets.
Bisection is particularly useful when interpolation requests are uncorrelated—meaning there is no pattern in the sequence of target points that could be exploited for faster searching.

```{code-cell} ipython3
def bisection(xs, target):
    l, h = 0, len(xs) - 1
    while h - l > 1:
        m = (h + l) // 2
        if target >= xs[m]:
            l = m
        else:
            h = m
    return l # returns index of the closest value less than or equal to target
```

The above function efficiently narrows down the interval to locate the index of the nearest value.

We can perform some tests:

```{code-cell} ipython3
for _ in range(10):
    xs = np.sort(np.random.uniform(0, 100, 10))
    v  = np.random.uniform(min(xs), max(xs))
    i  = bisection(xs, v)
    print(f'{xs[i]} <= {v} < {xs[i+1]}')
```

### Hunting Method

For cases where interpolation points are close together in sequence---common in applications with gradually changing target values---the hunting method offers faster performance than bisection from scratch.
Hunting takes advantage of the idea that, if the previous interpolation point is nearby, the search can start close to the last found position and "hunt" outward in expanding steps to bracket the target value.
Once the bracket is located, the search is refined using a quick bisection.

The hunting method is beneficial for correlated data requests, where successive target values are close, as it can skip large portions of the data and converge faster than starting from scratch each time.

```{code-cell} ipython3
def hunt(xs, target, i_last):
    n = len(xs)
    assert 0 <= i_last < n - 1

    # Determine the search direction based on the target value
    if target >= xs[i_last]:
        l, h, step = i_last, min(n-1, i_last+1), 1
        while h < n - 1 and target > xs[h]:
            l, h = h, min(n-1, h+step)
            step *= 2
    else:
        l, h, step = max(0, i_last-1), i_last, 1
        while l > 0 and target < xs[l]:
            l, h = max(0, l-step), l
            step *= 2

    # Refine with bisection within the bracketed range
    return bisection(xs[l:h+1], target) + l
```

```{code-cell} ipython3
i = 5
for _ in range(10):
    xs = np.sort(np.random.uniform(0, 100, 10))
    v  = np.random.uniform(min(xs), max(xs))
    i  = hunt(xs, v, i)
    print(f'{xs[i]} <= {v} < {xs[i+1]}')
```

### Linear Interpolation Using the Hunting Method

Once the nearest position is identified, interpolation proceeds with the closest data points.
Here, we implement a simple linear interpolation using the hunting method to locate the starting position, then use it to calculate the interpolated value.

```{code-cell} ipython3
class Interpolator:
    def __init__(self, xs, ys):
        assert len(xs) == len(ys)
        self.xs, self.ys = xs, ys
        self.i_last = len(xs)//2

    def __call__(self, target, search_method='hunt'):
        if search_method == 'hunt':
            i = hunt(self.xs, target, self.i_last)
        elif search_method == 'bisection':
            i = bisection(self.xs, target)
        else:
            i = linear(self.xs, target)
        self.i_last = i  # Update last position for future hunts
            
        # Linear interpolation using the two nearest points
        x0, x1 = self.xs[i], self.xs[i + 1]
        y0, y1 = self.ys[i], self.ys[i + 1]

        return (y1 - y0) * (target - x0) / (x1 - x0) + y0
```

```{code-cell} ipython3
def f(x):
    return np.exp(-0.5 * x * x)
    
Xs = np.sort(np.random.uniform(-5, 5, 20))
Ys = f(Xs)

fi = Interpolator(Xs, Ys)
```

```{code-cell} ipython3
from matplotlib import pyplot as plt

xs = np.linspace(min(Xs), max(Xs), 100)
ys = np.array([fi(x) for x in xs])

plt.plot(xs, ys, '.-')
plt.plot(Xs, Ys, 'o')
```

Let's test if our claim in terms of performance works in real life.

```{code-cell} ipython3
Xs = np.sort(np.random.uniform(-5, 5, 100))
Ys = f(Xs)
fi = Interpolator(Xs, Ys)

xs = np.linspace(min(Xs), max(Xs), 10000)
```

```{code-cell} ipython3
%timeit ys = np.array([fi(x, search_method='linear') for x in xs])
```

```{code-cell} ipython3
%timeit ys = np.array([fi(x, search_method='bisection') for x in xs])
```

```{code-cell} ipython3
%timeit ys = np.array([fi(x, search_method='hunt') for x in xs])
```

## Polynomial Interpolation and Extrapolation

Given $M$ data points $(x_0, y_0), (x_1, y_1), \dots, (x_{M-1}, y_{M_1})$, there exists a unique polynomial of degree $M-1$ that pass through all $M$ points exactly.

+++

### Lagrange's formula

This polynomial is given by Lagrange's classical formula,
\begin{align}
P_{M-1}(x)
&= \frac{(x-x_1)(x-x_2)\dots(x-x_{M-1})}{(x_0-x_1)(x_0-x_2)\dots(x_0-x_{M-1})} y_0 \\
&+ \frac{(x-x_0)(x-x_2)\dots(x-x_{M-1})}{(x_1-x_0)(x_1-x_2)\dots(x_1-x_{M-1})} y_1 + \dots \\
&+ \frac{(x-x_0)(x-x_2)\dots(x-x_{M-2})}{(x_{M-1}-x_0)(x_{M-1}-x_1)\dots(x_{M-1}-x_{M-2})} y_{M-1}
\end{align}
Using summation and product notations, one may rewrite Lagrange's formula as
\begin{align}
P_{M-1}(x)
= \sum_{m=0}^{M-1} \frac{\prod_{n=0,n\ne m}^{M-1}(x-x_n)}{\prod_{n=0,n\ne m}^{M-1}(x_m-x_n)} y_m
\end{align}
Substituting $x = x_{m'}$ for $0 \le m; < M$, it is straightforward to show
\begin{align}
P_{M-1}(x_{m'})
= \sum_{m=0}^{M-1} \delta_{mm'} y_m
\end{align}
and hence $P_{M-1}(x)$ does pass through all data points.

+++

### Neville's Algorithm

Although one may directly implement Lagrange's formula, it does not offer a way to estimate errors.
Instead, we will use Neville's algorithm, which constructs an interpolating polynomial by combining values in a recursive manner.
This approach avoids some of the issues in Lagrange interpolation and is particularly useful for small sets of points where we need an error estimate along with the interpolation result.

Although Numerical Receipts usually give excellent explainations of numerical methods, its section on Neville's Algorithm is a bit confusing.
Here, we try to use some python codes to motivate the algorithm step by step.

+++

1. Note that a polynomial of 0 dgree is simply a constant.
   We use $P_m$ to denote the 0 degree polynomails that approxmation points $(x_m, y_m)$.
   Hence, $P_m = y_m$.
   This is represented by the horizontal bars in the following figure.

```{code-cell} ipython3
Xs = np.sort(np.random.uniform(-5, 5, 100))
Ys = f(Xs)

plt.scatter(Xs, Ys, marker='_', label=r'$P_m$: polynomials with 0 degree')
plt.xlim(-0.5, 0.5)
plt.ylim( 0.9, 1.05)
plt.legend()
```

2. To improve the accuracy of the approximation, we try to linearly interpolate two nearby points $(x_{m'}, y_{m'})$ and $(x_{m'+1}, y_{m'+1})$.
   For book keeping reason, we will call this polynomial of 1 degree $P_{m',m'+1}$.
   Recall the previous definition $P_{m'} = y_{m'}$ and $P_{m'+1} = y_{m'+1}$,
   we may now use the "two-point form" and write down:
   \begin{align}
   \frac{P_{m',m'+1} - P_{m'}}{x - x_{m'}} &= \frac{P_{m'+1} - P_{m'}}{x_{m'+1} - x_{m'}} \\
   P_{m',m'+1} - P_{m'} &= \frac{x - x_{m'}}{x_{m'+1} - x_{m'}}(P_{m'+1} - P_{m'}) \\
   P_{m',m'+1} &= P_{m'} + \frac{x - x_{m'}}{x_{m'+1} - x_{m'}}(P_{m'+1} - P_{m'}) \\
   P_{m',m'+1} &= \frac{x_{m'+1} - x_{m'}}{x_{m'+1} - x_{m'}}P_{m'} + \frac{x - x_{m'}}{x_{m'+1} - x_{m'}}(P_{m'+1} - P_{m'}) \\
   P_{m',m'+1} &= \left(\frac{x_{m'+1} - x_{m'}}{x_{m'+1} - x_{m'}} - \frac{x - x_{m'}}{x_{m'+1} - x_{m'}}\right)P_{m'} + \frac{x - x_{m'}}{x_{m'+1} - x_{m'}}P_{m'+1} \\
   P_{m',m'+1} &= \frac{x_{m'+1} - x}{x_{m'+1} - x_{m'}}P_{m'} + \frac{x - x_{m'}}{x_{m'+1} - x_{m'}}P_{m'+1} \\
   P_{m',m'+1} &= \frac{(x - x_{m'+1})P_{m'} + (x_{m'} - x)P_{m'+1} }{x_{m'} - x_{m'+1}}
   \end{align}
   This is nothing but a special case of equation (3.2.3) in Numerical Recipes 3rd Edition in C++.

```{code-cell} ipython3
Pmm1s = []
for m in range(len(Xs)-1):
    Pmm1s.append(lambda x: ((x - Xs[m+1]) * Ys[m] + (Xs[m] - x) * Ys[m+1]) / (Xs[m] - Xs[m+1]))
```

```{code-cell} ipython3
plt.scatter(Xs, Ys, marker='_', label=r'$P_m$: polynomials with 0 degree')

for m, Pmm1 in enumerate(Pmm1s):
    xs = np.linspace(Xs[m], Xs[m+1], 100)
    ys = Pmm1(xs)
    plt.plot(xs, ys)

plt.xlim(-0.5, 0.5)
plt.ylim( 0.9, 1.05)
plt.legend()
```

```{code-cell} ipython3
Pmm1m2s = []
for m in range(len(Xs)-2):
    Pmm1  = lambda x: ((x - Xs[m+1]) * Ys[m  ] + (Xs[m  ] - x) * Ys[m+1]) / (Xs[m  ] - Xs[m+1])
    Pm1m2 = lambda x: ((x - Xs[m+2]) * Ys[m+1] + (Xs[m+1] - x) * Ys[m+2]) / (Xs[m+1] - Xs[m+2])
    Pmm1m2s.append(
        lambda x: ((x - Xs[m+2]) * Pmm1(x) + (Xs[m] - x) * Pm1m2(x)) / (Xs[m] - Xs[m+2])
    )
```

```{code-cell} ipython3
plt.scatter(Xs, Ys, marker='o', label=r'$P_m$: polynomials with 0 degree')

for m, Pmm1 in enumerate(Pmm1s[:-1]):
    xs = np.linspace(Xs[m], Xs[m+2], 100)
    ys = Pmm1(xs)
    plt.plot(xs, ys)

for m, Pmm1m2 in enumerate(Pmm1m2s):
    xs = np.linspace(Xs[m], Xs[m+2], 100)
    ys = Pmm1m2(xs)
    plt.plot(xs, ys, ':')

plt.xlim(-0.5, 0.5)
plt.ylim( 0.9, 1.05)
plt.legend()
```

3. By the same token, to improve the accuracy of the approximation, we linearly interpolate $P_{m'',m''+1}$ and $P_{m''+1,m''+2}$.
   We will call this polynomial of 2 degrees $P_{m'',m''+1,m''+2}$:
   \begin{align}
   P_{m'',m''+1,m''+2} &= \frac{(x - x_{m''+2})P_{m'',m''+1} + (x_{m''} - x)P_{m''+1,m''+2} }{x_{m'} - x_{m'+2}}.
   \end{align}

+++

4. Doing this recursively, we obtain Neville's algorithm, equation (3.2.3) in Numerical Recipes:
   \begin{align}
   P_{m,m+1,\dots,m+n} &= \frac{(x - x_{m+n})P_{m,m+1,\dots,m+n-1} + (x_{m} - x)P_{m+1,m+2,\dots,m+n} }{x_{m} - x_{m+n}}.
   \end{align}

5. Recalling the **catastrophic cancellation** discussed in [](data.md), given that $P_{m,m+1,\dots,m+n} - P_{m,m+1,\dots,m+n-1}$ has the meaning of "small correction", it is better to keep track of this small quantities instead of $P_{m,m+1,\dots,m+n}$ themselves.
   Following Numerical Recipes, we define
   \begin{align}
   C_{n,m} &\equiv P_{m,m+1,\dots,m+n} - P_{m,m+1,\dots,m+n-1} \\
   D_{n,m} &\equiv P_{m,m+1,\dots,m+n} - P_{m+1,m+2,\dots,m+n}
   \end{align}
   Neville's algorithm can now be rewritten as
   \begin{align}
   D_{n+1,m} &= \frac{x_{m+n+1}-x}{x_m - x_{m+n+1}}(C_{n,m+1} - D_{n,m}) \\
   C_{n+1,m} &= \frac{x_{n}-x}{x_m - x_{m+n+1}}(C_{n,m+1} - D_{n,m})
   \end{align}
   From this expression, it is now clear that the $C$'s and $D$'s are the corrections that make the interpolation one order higher.

6. The final polynomial $P_{0,1,\dots,M-1}$ is equal to the sum of *any* $y_i$ plus a set of $C$'s and/or $D$'s that form a path through the family tree of $P_{m,m+1,\dots,m+n}$.

```{code-cell} ipython3
class PolynomialInterpolator:
    def __init__(self, xs, ys, n=None):
        if n is None:
            n = len(xs)
        
        assert len(xs) == len(ys)        
        assert len(xs) >= n
        
        self.xs, self.ys, self.n = xs, ys, n

    def __call__(self, target, search_method='hunt'):

        C = np.copy(self.ys)
        D = np.copy(self.ys)
        
        i = np.argmin(abs(self.xs - target))
        y = self.ys[i]
        i-= 1
 
        for n in range(1,self.n):
            ho  = self.xs[:-n] - target
            hp  = self.xs[+n:] - target
            w   = C[1:self.n-n+1] - D[:-n]
            den = ho - hp
            if any(den == 0):
                raise Exception("two input xs are (to within roundoﬀ) identical.")
            else:
                f = w / den
            D[:-n] = hp * f
            C[:-n] = ho * f
                
            if 2*(i+1) < (self.n-n):
                self.dy = C[i+1]
            else:
                self.dy = D[i]
                i -= 1

            y += self.dy
        
        return y
```

```{code-cell} ipython3
Xs = np.linspace(0,2*np.pi,10)
Ys = np.sin(Xs)
P = PolynomialInterpolator(Xs, Ys)

xs = np.linspace(0,2*np.pi,315)
ys = []
es = []
for x in xs:
    ys.append(P(x))
    es.append(P.dy)
ys = np.array(ys)
es = np.array(es)

fig, axes = plt.subplots(2,1,figsize=(8,6))
axes[0].scatter(Xs, Ys)
axes[0].plot(xs, ys, '-', color='r')
axes[1].semilogy(xs, abs(es))
```

What will happen if we try to extrapolate?

+++

## Cubic Spline Interpolation

## Rational Function Interpolation and Extrapolation

## Coefficients of the Interpolating Polynomial

## Interpolation on a Grid in Multidimensions

## Interpolation on Scattered Data in Multidimensions

## Laplace Interpolation

## Conclusion and Discussion
