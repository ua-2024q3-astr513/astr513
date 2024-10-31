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

### Newton-Raphson Method for Nonlinear Systems of Equations

The Newton-Raphson Method is also a powerful technique for solving systems of nonlinear equations, where variables are interdependent in complex ways.
In a system of linear equations, we can directly solve for variables by inverting the coefficient matrix.
However, for nonlinear systems, such direct inversion is not possible.
Instead, the Newton-Raphson method provides an iterative approach, refining an initial guess using the Jacobian matrix of partial derivatives to approximate solutions.

In astrophysics, this method is essential for modeling complex phenomena such as stellar structure, orbit determination, and radiative transfer.
For example, in stellar structure modeling, equations governing temperature, pressure, and density must be solved simultaneously to describe a star's internal equilibrium.
Similarly, in orbit determination, nonlinear gravitational interactions are iteratively solved to yield accurate trajectories of celestial bodies.

+++

The Newton-Raphson method for systems proceeds as follows:
1.  Define Initial Guess: Start with an initial vector $\mathbf{x}^{(0)}$.
2.  Compute the Jacobian Matrix $J$:
    This matrix of partial derivatives represents the sensitivity of each equation with respect to each variable:
    \begin{align}
    J(\mathbf{x}) = \begin{bmatrix}
    \frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \dots & \frac{\partial f_1}{\partial x_n} \\
    \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \dots & \frac{\partial f_2}{\partial x_n} \\
    \vdots & \vdots & \ddots & \vdots \\
    \frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \dots & \frac{\partial f_n}{\partial x_n}
    \end{bmatrix}
    \end{align}
3.  Update Solution: Solve the linear system $J(\mathbf{x}) \Delta\mathbf{x} = -\mathbf{F}(\mathbf{x})$ for the correction term $\Delta \mathbf{x}$ and update $\mathbf{x}$ as:
    \begin{align}
    \mathbf{x}^{(n+1)} = \mathbf{x}^{(n)} + \Delta \mathbf{x}
    \end{align}
4.  Check Convergence: Iterate until $\|\Delta \mathbf{x}\|$ or $\|\mathbf{F}(\mathbf{x})\|$ is below a specified tolerance.

```{code-cell} ipython3
import numpy as np

# Newton-Raphson method for systems
def newton_system(F, J, X0, tol=1e-6, max_iter=100):
    for _ in range(max_iter):
        F0 = F(X0)
        J0 = J(X0)
        dX = np.linalg.solve(J0, -F0)
        X  = X0 + dX
        if np.linalg.norm(dX) < tol:
            return X
        X0 = X
    raise ValueError("Maximum iterations reached without convergence")
```

Consider a simplified model with two hypothetical equations that could represent balance conditions in astrophysics:
\begin{align}
\begin{cases}
f_1(x, y) = x^2 + y^2 - 4 = 0 \\
f_2(x, y) = e^x + y - 1 = 0
\end{cases}
\end{align}

```{code-cell} ipython3
# Define the function vector
def F(X):
    return np.array([X[0]**2 + X[1]**2 - 4, np.exp(X[0]) + X[1] - 1])

# Define the Jacobian matrix
def J(X):
    return np.array([[2 * X[0], 2 * X[1]], [np.exp(X[0]), 1]])

# Example usage
X0   = np.array([1.0, 1.0])
Root = newton_system(F, J, X0)

print("Approximate root:")
print("  R  = ",   Root )
print("F(R) = ", F(Root))
```

Similarly, we may take advantage of `JAX`'s autodifff to obtain the Jacobian:

```{code-cell} ipython3
from jax import jacfwd, numpy as jnp

# Newton-Raphson method for systems
def autonewton_system(F, X0, tol=1e-6, max_iter=100):
    J = jacfwd(F)
    for _ in range(max_iter):
        F0 = F(X0)
        J0 = J(X0)
        dX = np.linalg.solve(J0, -F0)
        X  = X0 + dX
        if np.linalg.norm(dX) < tol:
            return X
        X0 = X
    raise ValueError("Maximum iterations reached without convergence")
```

```{code-cell} ipython3
# Define the function vector with `jnp`
def F(X):
    return jnp.array([X[0]**2 + X[1]**2 - 4, jnp.exp(X[0]) + X[1] - 1])

# Example usage
X0   = np.array([1.0, 1.0])
Root = autonewton_system(F, X0)

print("Approximate root:")
print("  R  = ",   Root )
print("F(R) = ", F(Root))
```

This section introduced multiple root finding methods.
We covered the Bisection Method, a reliable but slower approach, and the more efficient Newton-Raphson Method, which converges quadratically but requires the function's derivative.
We then discuss the Secant Method, which approximates the derivative, allowing it to bypass the need for explicit derivative calculations.
For complex applications, modern libraries like `JAX` provide tools to automatically differentiate functions, making methods like Newton-Raphson easier to implement even for complex functions.
These techniques are widely applicable in astrophysics, where precise root approximations are needed.

+++

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
