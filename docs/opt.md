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

Optimization methods aim to find values that maximize or minimize an objective function, making them essential across fields like engineering, economics, and data science.
In fact, the whole action principle is about optimization.

In astronomy and astrophysics, optimization is crucial for model fitting, data analysis, and instrument calibration.
For example, in astrophysical model fitting, researchers use optimization to adjust parameters in models of, e.g., AGNs or galaxies, minimizing the difference between observed data and theoretical predictions.
This allows them to infer physical properties like energy output and elemental composition.

Optimization is also essential in orbit determination and mission planning, where it helps calculate fuel-efficient trajectories for spacecraft navigating complex gravitational fields to reach distant targets.
Precise trajectory optimization ensures mission success while conserving fuel resources.

Optimization underpins machine learning applications in astronomy as well, where models are trained by minimizing a loss function that measures prediction errors.
It enables tasks like galaxy classification and transient detection across large astronomical datasets.

Together, these applications demonstrate how optimization facilitates efficient, high-precision solutions in traditional fields and modern astronomy, advancing data quality, model accuracy, and scientific discovery.

+++

### Gradient Descent Methods

Gradient Descent is one of the most widely used optimization techniques, particularly effective for high-dimensional problems in fields like machine learning and astrophysics.
The method iteratively moves toward the minimum of a function by taking steps proportional to the negative of the gradient (or slope) at the current point, guiding the process toward the lowest value of the function.
For differentiable objective functions, gradient descent is fundamental in minimizing errors, making it essential for training machine learning models.

For a function $f(x)$, the gradient $\nabla f(x)$ points in the direction of steepest ascent, so taking steps in the opposite direction (the negative gradient) reduces the function's value. Gradient descent iteratively updates the parameters $x$ as follows:
\begin{align}
x_{n+1} = x_n - \alpha \nabla f(x_n)
\end{align}
where $\alpha$ is the learning rate, which controls the step size.
The choice of $\alpha$ is crucial to convergence: a large $\alpha$ may cause divergence (where the steps overshoot the minimum and fail to settle), while a very small $\alpha$ leads to slow convergence.
Proper tuning of $\alpha$ is therefore essential to ensure the algorithm successfully converges to a minimum without oscillating or diverging.

```{code-cell} ipython3
def gradient_descent(df, x, alpha, imax):
    for _ in range(imax):
        x -= alpha * df(x)
    return x
```

```{code-cell} ipython3
# Define the function and its gradient
def f(x):
    return (x - 3)**2 + 4

def df(x):
    return 2 * (x - 3)

# Parameters for gradient descent
x0    = 0.0  # Starting point for optimization
alpha = 0.1
imax  = 100

# Run gradient descent
xmin = gradient_descent(df, x0, alpha, imax)
print("Approximate minimum:")
print("  xmin  = ",   xmin )
print("f(xmin) = ", f(xmin))
```

We may visualize it by keep track of the history.
What will happen if we choose a large learning rate?

```{code-cell} ipython3
def gradient_descent_hist(df, x, alpha, imax):
    X = [x]
    for _ in range(imax):
        X.append(X[-1] - alpha * df(X[-1]))
    return X
```

```{code-cell} ipython3
from matplotlib import pyplot as plt

X = np.linspace(0, 6, 6001)
plt.plot(X, f(X))

alpha = 0.9

X = np.array(gradient_descent_hist(df, x0, alpha, imax))
print(X[-1])

plt.plot(X, f(X), '-o')
plt.xlim(2.5, 3.5)
plt.ylim(3.95,4.3)
```

Similar to our implementation of Newton-Raphson Method, it is possible to employ `JAX` to automatically obtain the derivative.
Here is an updated version of automatic gradient descent.

```{code-cell} ipython3
def autogd_hist(f, x, alpha, imax):
    df = grad(f)
    X  = [x]
    for _ in range(imax):
        X.append(X[-1] - alpha * df(X[-1]))
    return X
```

```{code-cell} ipython3
# Define the function and its gradient
def f(x):
    return (x - 3)**2 + 4

# Parameters for gradient descent
x0    = 0.0  # Starting point for optimization
alpha = 0.9
imax  = 100

# Run gradient descent
Xmin = np.array(autogd_hist(f, x0, alpha, imax))
print("Approximate minimum:")
print("  xmin  = ",   Xmin[-1] )
print("f(xmin) = ", f(Xmin[-1]))

X = np.linspace(0, 6, 6001)
plt.plot(X,    f(X))
plt.plot(Xmin, f(Xmin), '-o')
plt.xlim(2.5, 3.5)
plt.ylim(3.95,4.3)
```

### Gradient Descent with `JAX` for Multiple Dimensions

Multidimensional gradient descent is crucial for optimizing functions involving multiple parameters, making it essential for applications like model fitting and deep learning.
In model fitting tasks across scientific domains, such as astrophysics, gradient descent helps refine models by iteratively adjusting parameters to minimize discrepancies between observed data and theoretical predictions.
For example, in galaxy modeling, each parameter in the function may represent a physical property, such as brightness, size, or position, and gradient descent enables efficient optimization of these parameters to achieve the best fit to observational data.

In deep learning, multidimensional gradient descent is foundational, as modern neural networks can contain millions of parameters.
During training, gradient descent is used to minimize a loss function that measures the error between the model's predictions and actual outcomes.
Automatic differentiation with `JAX` simplifies the calculation of these gradients, allowing deep learning practitioners to easily train complex models without manually computing derivatives.
This automatic handling of gradient calculations is particularly valuable when networks have intricate architectures, such as convolutional or recurrent layers, where gradients must be calculated for a vast number of parameters across multiple dimensions.

The following example demonstrates how to use `JAX` to perform gradient descent on a multivariable function and track the optimization path, illustrating how gradient descent converges toward a minimum.

+++

In this example, we will minimize the function:
\begin{align}
f(x, y) = (x - 3)^2 + (y + 4)^2
\end{align}
where the minimum is at $(x, y) = (3, -4)$.
By tracking each update step, we can visualize the optimization path as it approaches the minimum.

```{code-cell} ipython3
# Function to perform gradient descent with history tracking
def autogd_hist(f, X, alpha, imax):
    df = grad(f)  # Use JAX to compute gradient
    Xs = [np.array(X)]
    for _ in range(imax):
        Xs.append(Xs[-1] - alpha * df(Xs[-1]))  # Gradient descent update
    return jnp.array(Xs)
```

```{code-cell} ipython3
# Define a multivariable function
def f(X):
    x, y = X
    return (x - 3)**2 + 2 * (y + 4)**2

# Parameters for gradient descent
X0    = jnp.array([0.0, 0.0]) # Starting point for optimization
alpha = 0.1                   # Learning rate
imax  = 100                   # Number of iterations

# Run gradient descent with history tracking
Xs = autogd_hist(f, X0, alpha, imax)
print("Approximate minimum:")
print("  xmin  =",   Xs[-1] )
print("f(xmin) =", f(Xs[-1]))

# Plot the function and gradient descent path
x_vals = jnp.linspace(-1, 7, 100)
y_vals = jnp.linspace(-8, 0, 100)
X, Y   = jnp.meshgrid(x_vals, y_vals)
Z      = f([X, Y])

plt.contour(X, Y, Z, levels=20)
plt.plot(Xs[:,0], Xs[:,1], '-o', color='red')
plt.xlabel('x')
plt.ylabel('y')
plt.gca().set_aspect('equal')
```

To demonstrate a more complex optimization scenario, let's consider fitting a multi-parameter model to noisy data.
We will use polynomial regression as our example, where we fit a polynomial curve to data points by optimizing the coefficients.
This is a non-trivial problem because, as the degree of the polynomial increases, the number of parameters grows, resulting in a high-dimensional optimization task.

```{code-cell} ipython3
groundtruth = np.array([1.2, -3, 0.5, 1.0, -1.8, 2.0, -0.1])

Xdata = np.linspace(-1, 1, 100)
Ytrue = sum(c * Xs**i for i, c in enumerate(groundtruth))
Ydata = Ytrue + np.random.normal(scale=0.1, size=Xdata.shape)
```

```{code-cell} ipython3
plt.plot(Xdata, Ytrue)
plt.plot(Xdata, Ydata)
```

```{code-cell} ipython3
# Define polynomial model
def model(Xs, Cs):
    return sum(c * Xs**i for i, c in enumerate(Cs))
    
# Define the objective function
def chi2(Cs):
    Ymodel = model(Xdata, Cs)
    return jnp.mean((Ymodel - Ydata)**2)

# Parameters for gradient descent
C0    = jnp.zeros(len(groundtruth)) # Start with zeros as initial coefficients
alpha = 0.1                         # Learning rate
imax  = 1000                        # Number of iterations

Cs = autogd_hist(chi2, C0, alpha, imax)

print("Optimized coefficients:", Cs[-1])
print("True coefficients:",      groundtruth)
print("Mean Squared Error:",     np.mean((groundtruth - Cs[-1])**2))
```

```{code-cell} ipython3
plt.scatter(Xdata, Ydata, color='blue', label='Noisy Data', alpha=0.5)
plt.plot(Xdata, Ytrue, 'g--', label='True Polynomial')
skip = 100
for i, Ci in enumerate(Cs[::skip]):
    Yfit = model(Xdata, Ci)
    plt.plot(Xdata, Yfit, 'r', alpha=skip*i/imax, label='Fitted Polynomial' if skip*i == imax else '')
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
```

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
