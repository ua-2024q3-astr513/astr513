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

# Integration of Ordinary Differential Equations. I

In this lecture, we will explore numerical methods for solving ordinary differential equations (ODEs), focusing on initial value problems (IVPs).
We will begin by discussing the importance of integrating differential equations in physics and the limitations of analytical solutions.
We will then introduce Euler's method, both forward (explicit) and backward (implicit) forms.
We will examine how truncation errors and round-off errors affect these methods.
Additionally, we will demonstrate that the backward Euler method is more stable than the forward Euler method by analyzing their stability regions.

+++

## Introduction to Numerical ODEs

Ordinary differential equations (ODEs) are foundational to understanding and predicting the behavior of physical systems.
Many of the dynamic processes that occur in nature, from the motion of planets under gravitational forces to the oscillations of atoms in a crystal lattice, are described by differential equations.
For instance, Newton's second law, which states that the force on an object is equal to its mass times its acceleration, can be written as:
\begin{align}
F = m a = m \frac{d^2x}{dt^2}
\end{align}
This equation, along with others describing physical laws, typically involves derivatives of a function, representing rates of change, which relate the state of a system at one moment to its state at another.
In cases where physical fields, like the electromagnetic or gravitational fields, are described, partial differential equations (PDEs) are used since they involve rates of change with respect to multiple variables.
However, in this lecture, we focus on ODEs where functions depend on a single variable, typically time.

It is not always possible to solve ODEs analytically
Analytic solutions are often limited to simple or highly idealized cases.
Real-world systems tend to involve complex, nonlinear relationships that cannot be resolved through straightforward integration or algebraic manipulation.
For these cases, numerical methods provide a powerful alternative.

+++

### Problem Definition and Types of Differential Equations

Consider the following two types of differential equations, which differ in their dependence on variables:

1.  The first type has the form:
    \begin{align}
    \frac{dx}{dt} = f(t)
    \end{align}
    Here, the right-hand side (RHS) of the equation, $f(t)$, depends only on the independent variable $t$ and not on $x$, the dependent variable we are solving for.
2.	The second type is given by:
    \begin{align}
    \frac{dx}{dt} = f(x, t)
    \end{align}
    In this equation, the RHS depends on both $t$ and $x$, making it more complicated since $x$, the quantity we are trying to determine, appears in the expression that defines its rate of change.

+++

The first type of differential equation, where $f(t)$ depends only on $t$, can be solved by direct integration.
This is because we can integrate $f(t)$ with respect to $t$ to find $x(t)$:
\begin{align}
x(t) = x(t_0) + \int_{t_0}^t f(t') \, dt'
\end{align}
This approach is feasible when $f(t)$ is a function that can be integrated analytically.
In cases where this is not possible, numerical integration techniques, such as the trapezoidal rule or Simpson's rule, can be used to approximate the solution.
We covered these numerical integration techniques in an earlier lecture on "[Numerical Integration of Functions](integration.md)".

+++

The second type of differential equation, where $f(x, t)$ depends on both $x$ and $t$, is more complicated.
We cannot solve this type by direct integration because $x$, the function we are trying to find, is also part of the expression on the RHS. For example, attempting to directly integrate is not feasible since $x$ is unknown within the integral itself.
In other words, we do not know $x$ at intermediate points between $t_0$ and $t$, so we cannot compute the integral without first determining $x(t)$ at these points.

This is where numerical methods for ODEs become essential.
Instead of solving for $x(t)$ in one go, we use numerical methods to approximate the solution by advancing in small increments over the interval of interest.
This way, we iteratively approximate $x(t)$ at discrete points in time, allowing us to handle equations where direct integration is not possible.

Since many real-world systems fall into the second category of differential equations, where the dependence on  $x$ and $t$ is nonlinear, numerical methods are widely used in scientific computing.
These methods break down the problem into manageable steps and enable us to approximate solutions even when analytic solutions do not exist or are difficult to obtain.
Numerical methods not only provide practical solutions to ODEs but also allow for flexibility in modeling a wide range of phenomena, such as chaotic systems, non-linear oscillations, and biological systems with interacting populations.

+++

## Euler's Method

Euler's method is the simplest techniques for numerically solving ordinary differential equations.
This method provides an easy way to approximate the solution of an initial value problem (IVP) by advancing one small step at a time.
We can apply Euler's method to an ODE of the form:
\begin{align}
\frac{dx}{dt} = f(x, t), \quad x(t_0) = x_0
\end{align}
where $x_0$ is the initial value of $x$ at time $t = t_0$.
However, as we will see below, it is usually not recommanded in pratical calculations because of its stability properties.

+++

### Forward (Explicit) Euler Method

There are three simple ways to derive Euler's method.

The easiest way is simply hold $x$ fixed in $f(x, t)$ and apply the left Reimann sum.
The left Reimann sum is first order in step size by approximating $f(x, t)$ as a constant.
In this sense, holding $x$ is somewhat "self-consistent" in terms of approximation order.

Second, we recall the definition of a deriviative:
\begin{align}
  f(x, t) = \frac{dx}{dt} = \lim_{\Delta t\rightarrow 0}\frac{x(t + \Delta t) - x(t)}{\Delta h}.
\end{align}
If we simply remove the limit and keep the "finite difference" part, then it is trivial to show
\begin{align}
  x(t + \Delta t) &\approx x(t) + f(x(t), t)\Delta t.
\end{align}
Which is nothing but again the forward Euler method.
While very intuitive, the above two derivations do not formally show the order of the Euler method.

We may also consider a numerical approximation to the solution of an ODE.
We approximate the solution at time $t_{n+1} = t_n + \Delta t$ by using the Taylor expansion:
\begin{align}
x(t_{n+1}) = x(t_n) + f(x(t_n), t_n) \Delta t + \mathcal{O}(\Delta t^2)
\end{align}
Neglecting the higher-order terms in the expansion, we obtain once again the Forward Euler Method:
\begin{align}
x_{n+1} = x_n + f(x_n, t_n) \Delta t
\end{align}
It is thus a step-by-step approach that proceeds by evaluating $f(x, t)$ at each time point and then advancing to the next point.
It is an explicit method in 1st order.

+++

### Example

Let's consider solving the simpmle differential equation:
\begin{align}
\frac{dx}{dt} = \lambda x(t)
\end{align}

This equation has solution
\begin{align}
x(t) = \exp[\lambda(t-t_0)]
\end{align}

If we choose $\lambda = 1$ and $x(0) = 1$, the solutoin reduces to $x(t) = \exp(t)$.

```{code-cell} ipython3
# Let's visualize the solution:

import numpy as np
from matplotlib import pyplot as plt

Tp = np.linspace(0, 2, 2001)
Xp = np.exp(Tp)

plt.plot(Tp, Xp)
```

```{code-cell} ipython3
# Let's implement Euler's method, with history

def forwardEuler(f, x, t, dt, n):
    T = [t]
    X = [x]
    for i in range(n):
        t += dt
        x += f(x) * dt
        T.append(t)
        X.append(x)
    return T, X
```

```{code-cell} ipython3
# Let's test Euler's method

def f(x):
    return x

T, X = forwardEuler(f, 1, 0, 0.1, 20)

plt.plot(Tp, Xp)
plt.plot(T,  X, 'o-')
plt.xlabel('t')
plt.ylabel('x')
```

```{code-cell} ipython3
# The accuracy should increase when we use a small step

T, X = forwardEuler(f, 1, 0, 0.01, 200)

plt.plot(Tp, Xp)
plt.plot(T,  X)
plt.xlabel('t')
plt.ylabel('x')
```

```{code-cell} ipython3
# As always, we can study the convergence of the numerical method

def error(N=200):
    T, X = forwardEuler(f, 1, 0, 2/N, N)
    Xp = np.exp(T)
    return np.max(abs(X - Xp))

N = np.array([64, 128, 256, 512, 1024])
E = np.array([error(n) for n in N])

plt.loglog(N, 16/N, label='1/N')
plt.loglog(N, E, 'o:', label='Forward Euler')
plt.xlabel('N')
plt.ylabel(r'$\text{err} = max|x_\text{numeric} - x|$')
plt.legend()
```

### System of ODEs

In computational astrophysics, we often encounter systems governed by Newton's laws:
\begin{align}
m \frac{d^2 x}{dt^2} = f(x, t)
\end{align}

This equation is a second-order ordinary differential equation because it involves the second derivative of $x$ with respect to $t$.
However, it is often more practical to convert second-order ODEs into a system of first-order ODEs.
To do this, we introduce a new variable, $v = dx/dt$, representing the velocity of the object.
This substitution allows us to rewrite the second-order ODE as two coupled first-order equations:
\begin{align} 
\frac{dx}{dt} &= v\\
\frac{dv}{dt} &= \frac{1}{m}f(x, t)
\end{align}

This formulation provides a convenient framework for numerical methods, which are generally well-suited to solving systems of first-order ODEs.
To further simplify, we can express these equations in vector notation by defining $Y = [x, v]^t$ and $F = [v, f/m]^t$. The system then becomes:
\begin{align}
\frac{dY}{dt} = F(Y, t).
\end{align}

This vector form emphasizes the structure of the system and enables us to apply general numerical techniques to solve both equations simultaneously.

+++

To illustrate this approach, let's consider a classic example: the simple pendulum under gravity.
The motion of a pendulum of length $l$, swinging under gravity $g$, can be described by the second-order equation:
\begin{align}
\frac{d^2\theta}{dt^2} + \frac{g}{l} \sin\theta = 0
\end{align}

Here, $\theta(t)$ is the angle of the pendulum with the vertical, and the term $\sin \theta$ introduces nonlinearity, which makes the equation challenging to solve analytically.
Converting this equation into a system of first-order ODEs allows us to handle it more effectively with numerical methods.
We define $\Omega = \frac{d\theta}{dt}$, the angular velocity, leading to the following system:
\begin{align}
\frac{d\theta(t)}{dt} &= \Omega(t)\\
\frac{d\Omega(t)}{dt} &= - \frac{g}{l}\sin\theta(t)
\end{align}

In vector notation, we represent the system as:
\begin{align}
\frac{dY(t)}{dt} = F(Y, t)
\end{align}
where 
\begin{align}
Y &= \begin{bmatrix} \theta(t) \\ \Omega(t) \end{bmatrix}\\
F(Y, t) &= \begin{bmatrix} \Omega(t) \\ -\frac{g}{l} \sin \theta(t) \end{bmatrix}
\end{align}

+++

In later part of the lecture, we will try to solve this problem.
But to derive and compare different numerical methods, let's first reduce the problem to something that has analytical solutions.
Specifically, we can simplify the pendulum problem further by assuming small oscillations, where the angle $\theta$ is small enough that $\sin \theta \approx \theta$.
This approximation linearizes the equation of motion, reducing the system to a simple harmonic oscillator.
In this approximation, the equation of motion becomes:
\begin{align}
\frac{d^2 \theta}{dt^2} + \frac{g}{l} \theta = 0
\end{align}
As a result, the system of ODEs becomes:
\begin{align}
Y &= \begin{bmatrix} \theta(t) \\ \Omega(t) \end{bmatrix}\\
F(Y, t) &= \begin{bmatrix} \Omega(t) \\ -\frac{g}{l} \theta(t) \end{bmatrix}
\end{align}

```{code-cell} ipython3
# Let's first plot the analytical solution

Tp     = np.linspace(0, 10, 1001)
Thetap = 0.01 * np.sin(Tp)

plt.plot(Tp, Thetap)
```

```{code-cell} ipython3
# Thanks to operator overriding,
# our forward Euler method is almost ready to solve system of ODEs

def forwardEuler(f, x, t, dt, n):
    T = np.array(t)
    X = np.array(x)
    for i in range(n):
        t += dt
        x += np.array(f(*x)) * dt
        T = np.append( T, t )
        X = np.vstack((X, x))
    return T, X
```

```{code-cell} ipython3
# Compare the analytical and numerical solutions

def f(theta, omega):
    return omega, -theta

T, X = forwardEuler(f, (0, 0.01), 0, 0.01, 1000)

Theta = X[:,0]
Omega = X[:,1]

plt.plot(Tp, Thetap)
plt.plot(T,  Theta)
```

```{code-cell} ipython3
# Again, we can study the convergence of the numerical method

def error(N=100):
    T, X = forwardEuler(f, (0, 0.01), 0, 10/N, N)
    Theta  = X[:,0]
    Thetap = 0.01 * np.sin(T)
    return np.max(abs(Theta - Thetap))

N = np.array([64, 128, 256, 512, 1024])
E = np.array([error(n) for n in N])

plt.loglog(N, 0.5/N, label='1/N')
plt.loglog(N, E, 'o:', label='Forward Euler')
plt.xlabel('N')
plt.ylabel(r'$\text{err} = max|x_\text{numeric} - x|$')
plt.legend()
```

### Can We Improve the Method?

The Forward Euler method is only first-order accurate, meaning its error decreases linearly with the step size $\Delta t$.
While it is simple to implement, this method's convergence rate is limited, and achieving higher accuracy requires very small steps, which can become computationally expensive.
This naturally raises the question: can we improve the convergence rate of our numerical solution, reducing the error more rapidly as we use smaller steps?
To explore this, we can draw inspiration from our previous work on numerical integration.
In that context, we observed that the midpoint (or central) Riemann sum converges faster than the left or right Riemann sums.
This suggests that a midpoint approach may also provide advantages in solving ODEs.

+++

### A Midpoint Method

One possible improvement is to propose a midpoint method that attempts to evaluate the function $f$ at the midpoint between steps.
Mathematically, this approach can be expressed as:
\begin{align}
X_{n+1} = X_n + f(X_{n+1/2}, t_{n+1/2}) \Delta t
\end{align}
However, a significant obstacle arises: the midpoint value $f_{n+1/2} \equiv f(X_{n+1/2}, t_{n+1/2})$ is unknown at step $n$.
We need to know the value of $X$ at the midpoint $t_{n+1/2}$ to use this method, but this value cannot be calculated without already knowing the future values of $X$.
This issue makes a straightforward midpoint method impractical for generic ODEs, where $f$ depends on both $X$ and $t$.

An exception occurs if $f$ depends only on $t$, as in $\frac{dx}{dt} = f(t)$;
in such cases, a true midpoint method is feasible, which is nothing but out middle Reimann sum.
However, for most ODEs, including those where $f$ depends on $X$, a different approach is necessary.

+++

### A Simple Idea: The Second-Order Runge-Kutta Method

To work around this issue, we can approximate the midpoint value instead of calculating it exactly.
It suffices to find an approximate solution for $X$ at the half-step, which we denote as $Y_{n+1/2}$.
One way to do this is to use the Forward Euler method to compute an estimated value at the midpoint.
Specifically, we can approximate $X$ at $t_{n+1/2}$ as:
\begin{align}
\tilde{X}_{n+1/2} = X_n + f(X_n, t_n) \frac{\Delta t}{2}
\end{align}
Using this half-step approximation, we then proceed with a full step to find $X_{n+1}$ by evaluating $f$ at the midpoint:
\begin{align}
X_{n+1} = X_n + f(\tilde{X}_{n+1/2}, t_{n+1/2}) \Delta t
\end{align}
This approach, known as the second-order Runge-Kutta method, provides a way to incorporate midpoint information, achieving second-order accuracy.
The second-order Runge-Kutta method improves convergence by leveraging an approximate midpoint, resulting in a more accurate solution than the first-order Euler method without requiring a prohibitively small step size.

```{code-cell} ipython3
def RK2(f, x, t, dt, n):
    T = np.array(t)
    X = np.array(x)
    
    for i in range(n):
        k1 = dt * np.array(f(*(x         )))
        k2 = dt * np.array(f(*(x + 0.5*k1)))
        
        t += dt
        x += k2
        
        T = np.append( T, t)
        X = np.vstack((X, x))
        
    return T, X
```

```{code-cell} ipython3
N=64

T, X = RK2(f, (0, 0.01), 0, 10/N, N)

plt.plot(T, 0.01*np.sin(T))
plt.plot(T, X[:,0], 'o')
```

```{code-cell} ipython3
def error2(N=100):
    T, X = RK2(f, (0, 0.01), 0, 10/N, N)
    Theta  = X[:,0]
    Thetap = 0.01 * np.sin(T)
    return np.max(abs(Theta - Thetap))

N = np.array([64, 128, 256, 512, 1024])
E2 = np.array([error2(n) for n in N])

print(E2[-1])

plt.loglog(N, 1/N**2,  label='1/N^2')
plt.loglog(N, E, 'o:', label='Forward Euler')
plt.loglog(N, E2,'o-', label='RK2')
plt.xlabel('N')
plt.ylabel(r'$\text{err} = max|x_\text{numeric} - x|$')
plt.legend()
```
