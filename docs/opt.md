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

The Bisection Method is guaranteed to converge to a root if the function is continuous on $[l, h]$ and $f(l)$ and $f(h)$ have different signs.
However, its convergence rate is relatively slow, decreasing the interval size by half with each iteration.
This results in a linear convergence rate (of error at fixed step).

+++

### Newton-Raphson Method

The Newton-Raphson Method is based on the concept of using the tangent line at a point on the curve of a function to approximate its root.
It leverages the Taylor series expansion to iteratively move closer to the root, achieving quadratic convergence when the initial guess is close to the actual root.

Consider a function $f(x)$ that we want to find the root.
Suppose we have an initial guess $x_0$ near the root.
We can expand $f(x)$ around this point $x_0$ using the Taylor series:
\begin{align}
f(x) = f(x_0) + f'(x_0)(x - x_0) + \frac{f''(x_0)}{2}(x - x_0)^2 + \cdots
\end{align}
For simplicity, we approximate the function linearly by ignoring the higher-order terms. This gives us a linear approximation:
\begin{align}
f(x) \approx f(x_0) + f'(x_0)(x - x_0)
\end{align}
We want to find the value of $x_0$ where $f(x_0) = 0$.
Therefore, 
\begin{align}
x \approx x_0 - \frac{f(x_0)}{f'(x_0)}.
\end{align}
We may call this next approximation $x_{n+1}$.
Thus, we define the iterative update as:
\begin{align}
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
\end{align}
This formula is the foundation of the Newton-Raphson Method.
It has a simple geometric interpretation.
Starting with an initial guess $x_0$, we compute the tangent line to the function $f(x)$ at $x_0$.
The tangent line provides a linear approximation of $f(x)$ near $x_0$, and where this tangent crosses the $x$-axis is our next estimate for the root, $x_1$.
By iterating this process, we continue to update our estimate, ideally converging quickly to the root.

+++

The Newton-Raphson Method is quadratically convergence near the root, which means that the error approximately squares with each iteration.
Specifically, if $x_n$ is close enough to the root, the error at the next step, $|x_{n+1} - x_\infty|$, is roughly proportional to $|x_n - x_\infty|^2$.
This quadratic convergence makes the Newton-Raphson Method highly efficient when close to the root.
However, it requires:
* **Non-Zero Derivative:** The method requires that $f'(x) \neq 0$ at each iteration.
  If $f'(x) = 0$ at any point, the update formula becomes undefined, and the algorithm will fail.
* **Good Initial Guess:** Convergence to the root is not guaranteed if the initial guess $x_0$ is too far from the actual root.
  Poor starting points can lead the method to diverge or to converge to the wrong root.
* **Well-Behaved Function:** The method performs best when $f(x)$ is smooth and continuous near the root.

+++

Here is a Python function that implements the Newton-Raphson Method with an initial guess and a tolerance level:

```{code-cell} ipython3
def newton(f, df, x0, tol=1e-6, imax=100):
    for _ in range(imax):
        f0, df0 = f(x0), df(x0)
        if df0 == 0:
            raise ValueError("Derivative is zero. No convergence.")
        x = x0 - f0 / df0
        if abs(x - x0) < tol:
            return x
        x0 = x
    raise ValueError("Maximum iterations reached without convergence")
```

```{code-cell} ipython3
f  = lambda x: x**3 - x - 2
df = lambda x: 3*x**2 - 1

initial_guess = 1.5
root = newton(f, df, initial_guess)
print("Approximate root:")
print("  x0  = ",   root )
print("f(x0) = ", f(root))
```

The Newton-Raphson Method is fast and efficient, especially when close to the root, due to its quadratic convergence.
However, in addition to the convergence conditions list above, it requires computing derivatives, making it less convenient compared to the bisection method.

Nevertheless, we learned about different derivative methods in this course, specifically, the machine accurate autodiff method.
By using it, we can implement a convenient Newton-Raphson Method:

```{code-cell} ipython3
import jax
jax.config.update("jax_enable_x64", True)

from jax import grad

def autonewton(f, x0, tol=1e-6, imax=100):
    df = grad(f)    
    for _ in range(imax):
        f0, df0 = f(x0), df(x0)
        if df0 == 0:
            raise ValueError("Derivative is zero. No convergence.")
        x = x0 - f0 / df0
        if abs(x - x0) < tol:
            return x
        x0 = x
    raise ValueError("Maximum iterations reached without convergence")
```

```{code-cell} ipython3
initial_guess = 1.5
root = autonewton(f, initial_guess)
print("Approximate root:")
print("  x0  = ",   root )
print("f(x0) = ", f(root))
```

### Secant Method

While in python we can use packages like `jax` to implement `autonewton()`, autodiff may not be available in some languages or in special embedded systems that has minimal compiler infrastructures.
In such a case, one may use the Secant Method, which approximates the derivative using values of $f(x)$ at two nearby points.
This makes it useful for functions where the derivative is difficult or expensive to compute.
The method is generally faster than the Bisection Method but can be less stable than Newton-Raphson if the initial points are not chosen carefully.

Here is a python implementation:

```{code-cell} ipython3
def secant(f, x0, x1, tol=1e-6, imax=100):
    for _ in range(imax):
        f0, f1 = f(x0), f(x1)
        if f0 == f1:
            raise ValueError("Division by zero in secant method.")
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        if abs(x2 - x1) < tol:
            return x2
        x0, x1 = x1, x2
    raise ValueError("Maximum iterations reached without convergence")
```

```{code-cell} ipython3
root = secant(f, 1, 2)
print("Approximate root:")
print("  x0  = ",   root )
print("f(x0) = ", f(root))
```

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
